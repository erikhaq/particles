#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include "osxbarrier.h"
#include <iostream>
#include <vector>

using namespace std;

//
//  global variables
//
int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;
CellMatrix cells;
pthread_mutex_t cs;
int numArrived = 0;

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

void update_cells_parallel(int start_row, int end_row, Particles &my_particles, CellMatrix &cells)
{
    for(int i = 0; i < my_particles.size(); i++)
    {
        particle_t *curr_particle = my_particles[i];
        Point p = get_cell_index(*curr_particle);
        if(p.y >= end_row || p.y < start_row)
        {
            pthread_mutex_lock(&cs);
            cells[p.y][p.x].push_back(curr_particle);
            pthread_mutex_unlock(&cs);
        }
        else
        {
            cells[p.y][p.x].push_back(curr_particle);
        }
    }
}

void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;
    int num_cells = cells.size();

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int rows_per_thread = (num_cells + n_threads - 1) / n_threads;

    int first_row = min(  thread_id     * rows_per_thread, num_cells);
    int last_row  = min( (thread_id+1)  * rows_per_thread, num_cells);

    int first = min(  thread_id    * particles_per_thread, n );
    int last  = min( (thread_id+1) * particles_per_thread, n );

    Particles my_particles;
    my_particles.clear();
    get_particles_from_rows(first_row, last_row, &my_particles, cells);

 
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        my_particles.clear();
        get_particles_from_rows(first_row, last_row, &my_particles, cells);
        
        //
        //  compute forces
        //
        for(int i = 0; i < my_particles.size(); i++)
        {
            particle_t *curr_particle = my_particles[i];
            curr_particle->ax = 0;
            curr_particle->ay = 0;
            apply_force(curr_particle, cells);            
        }

        
        pthread_barrier_wait( &barrier );
        
        //
        //  move particles
        //
        for(int i = 0; i < my_particles.size(); i++)
        {
            particle_t *curr_particle = my_particles[i];
            move(*my_particles[i]);            
        }
        
        pthread_barrier_wait( &barrier );
        clear_cells(first_row, last_row, cells);

        pthread_barrier_wait( &barrier );
        update_cells_parallel(first_row, last_row, my_particles, cells);
    
        pthread_barrier_wait( &barrier );
        
        //
        //  save if necessary
        //
        if( thread_id == 0 && fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    
    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  allocate resources
    //
    fsave = savename ? fopen( savename, "w" ) : NULL;

    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    int num_cells = get_num_cells();
    cells = CellMatrix(num_cells);
    init_cell_matrix(cells);
    update_cells(particles, cells, n);

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );

    pthread_mutex_init(&cs, NULL);
    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );
    
    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
