#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include "osxbarrier.h"
#include <iostream>

using namespace std;

//
//  global variables
//
int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;
CellMatrix cells;
pthread_mutex_t cellBarrier;
pthread_cond_t go;
int numArrived = 0;
//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}


/* a reusable counter barrier */
void cell_barrier() {
  pthread_mutex_lock(&cellBarrier);
  numArrived++;
  if (numArrived == n_threads) {
    numArrived = 0;
    update_cells(particles, cells, n);
    pthread_cond_broadcast(&go);
  } else
    pthread_cond_wait(&go, &cellBarrier);
  pthread_mutex_unlock(&cellBarrier);
}
//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;
    int num_cells = cells.size();

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int rows_per_thread = (num_cells + n_threads - 1) / n_threads;

    int first_row = min(  thread_id     * rows_per_thread, num_cells);
    int last_row  = min( (thread_id+1)  * rows_per_thread, num_cells);
    printf("first_row: %d\n", first_row);
    printf("last_row: %d\n", last_row);

    int first = min(  thread_id    * particles_per_thread, n );
    int last  = min( (thread_id+1) * particles_per_thread, n );
    // printf("first: %d\n", first);
    // printf("last: %d\n", last);
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        for( int i = first; i < last; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
        //         apply_force( particles[i], particles[j] );
            apply_force( &particles[i], cells);

        }
        
        pthread_barrier_wait( &barrier );
        
        //
        //  move particles
        //
        for( int i = first; i < last; i++ ) 
            move( particles[i] );
        
        // pthread_barrier_wait( &barrier );
        // clear_cells(first_row, last_row, cells);
        // printf("hej\n");


        // pthread_barrier_wait( &barrier );
        
        // printf("bajs\n");
        // update_cells_only(first, last, particles, cells);
        // pthread_barrier_wait( &barrier );
         cell_barrier();
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
    cout << num_cells << endl;
    cells = CellMatrix(num_cells);
    init_cell_matrix(cells);
    update_cells(particles, cells, n);


    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );
    pthread_mutex_init(&cellBarrier, NULL);
    pthread_cond_init(&go, NULL);
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
