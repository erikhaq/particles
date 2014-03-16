#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <omp.h>


void clear_row(int row, CellMatrix &cells)
{
    int num_cols = cells[row].size();
    for(int col = 0; col < num_cols; col++)
    {
        cells[row][col].clear();
    }

}
void update_cells_parallel(int start_row, int end_row, Particles &my_particles, CellMatrix &cells)
{
    for(int i = 0; i < my_particles.size(); i++)
    {
        Point p = get_cell_index(*my_particles[i]);
        if(p.y >= end_row || p.y < start_row)
        {
            #pragma omp critical
            cells[p.y][p.x].push_back(my_particles[i]);
            
        }
        else
        {
            cells[p.y][p.x].push_back(my_particles[i]);
        }
    }
}
//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    int num_cells = get_num_cells();
    CellMatrix cells(num_cells);
    init_cell_matrix(cells);
    update_cells(particles, cells, n);
    int n_threads = read_int( argc, argv, "-p", omp_get_max_threads() );
    omp_set_num_threads(n_threads);

    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int rows_per_thread = (num_cells + n_threads - 1) / n_threads;

        int first_row = min(  thread_id     * rows_per_thread, num_cells);
        int last_row  = min( (thread_id+1)  * rows_per_thread, num_cells);
        
        Particles my_particles;
        my_particles.clear();
        get_particles_from_rows(first_row, last_row, &my_particles, cells);

        for( int step = 0; step < 1000; step++ )
        {
            my_particles.clear();
            get_particles_from_rows(first_row, last_row, &my_particles, cells);
            
            //
            //  compute all forces
            //
            for(int i = 0; i < my_particles.size(); i++)
            {
                my_particles[i]->ax = 0;
                my_particles[i]->ay = 0;
                apply_force(my_particles[i], cells);
            }
            #pragma omp barrier            
            
            //
            //  move particles
            //
            for(int i = 0; i < my_particles.size(); i++)
            {
                move(*my_particles[i]);            
            }
            #pragma omp barrier

            clear_cells(first_row, last_row, cells);
            #pragma omp barrier

            update_cells_parallel(first_row, last_row, my_particles, cells);
            #pragma omp barrier
            
            //
            //  save if necessary
            //
            #pragma omp master
            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
    }
    
    simulation_time = read_timer( ) - simulation_time;
    
    printf("n = %d, \t threads = %d, \tsimulation time = %g seconds\n", n, n_threads,simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
