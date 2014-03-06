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
        particle_t *curr_particle = my_particles[i];
        Point p = get_cell_index(*curr_particle);
        if(p.y >= end_row || p.y < start_row)
        {
            #pragma omp critical
            cells[p.y][p.x].push_back(curr_particle);
            
        }
        else
        {
            cells[p.y][p.x].push_back(curr_particle);
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
    int n_threads = omp_get_max_threads();
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
        printf("first_row: %d\n", first_row);
        printf("last_row: %d\n", last_row);
        Particles my_particles;
        my_particles.clear();
        get_particles_from_rows(first_row, last_row, &my_particles, cells);
        printf("size: %d\n", my_particles.size());
        for( int step = 0; step < 1000; step++ )
        {
        //
        //  compute all forces
        //

            my_particles.clear();
            get_particles_from_rows(first_row, last_row, &my_particles, cells);
        
            //#pragma omp for
            for(int i = 0; i < my_particles.size(); i++)
            {
                particle_t *curr_particle = my_particles[i];
                curr_particle->ax = 0;
                curr_particle->ay = 0;
                apply_force( curr_particle, cells);            
            }
            #pragma omp barrier
            // for( int i = 0; i < n; i++ )
            // {
            //     particles[i].ax = particles[i].ay = 0;
            //     // for (int j = 0; j < n; j++ )
            //     //     apply_force( particles[i], particles[j] );
            //     apply_force(&particles[i], cells);
            // }
            
            //
            //  move particles
            //
            //#pragma omp for
            for(int i = 0; i < my_particles.size(); i++)
            {
                particle_t *curr_particle = my_particles[i];
                move(*my_particles[i]);            
            }
            #pragma omp barrier

            // for( int i = 0; i < n; i++ ) 
            //     move( particles[i] );
            
            // #pragma omp for
            // for(int r = 0; r < num_cells; r++)
            // {
            //     clear_row(r, cells);
            // }
            clear_cells(first_row, last_row, cells);
            #pragma omp barrier

            update_cells_parallel(first_row, last_row, my_particles, cells);
            #pragma omp barrier
            // #pragma omp for
            // for(int i = 0; i < n; i++)
            // {
            //     Point p = get_cell_index(particles[i]);
            //    // #pragma omp critical
            //     cells[p.y][p.x].push_back(&particles[i]);
            // }

            
            // #pragma omp single
            // update_cells(particles, cells, n);

            // #pragma omp single
            // update_cells_only(0, n, particles, cells);
            
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
