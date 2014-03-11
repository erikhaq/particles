#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <list>
#include "common.h"

using namespace std;

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    set_size(n);
    int num_cells = get_num_cells();
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    
    list<particle_t> top_reference;
    list<particle_t> bottom_reference;
    list<particle_t> my_particles;


    CellMatrix cells(num_cells);
   
    particle_t dummy;
    dummy.x = -1;
    dummy.y = -1;
    //
    //  set up the data partitioning across processors
    //
    
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    //
    //  allocate storage for local partition
    //
    

    int rows_per_thread = (num_cells + n_proc - 1) / n_proc;
    int top_row     = min(  rank     * rows_per_thread, num_cells); // the top row that the process is responsible for.
    int bottom_row  = min( (rank+1)  * rows_per_thread, num_cells); // the bottow row that this process needs but is not responsible for
    
    int my_amount;
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)

    list<particle_t> all_particles;
    all_particles.clear();
    init_cell_matrix(cells);
    all_particles.push_back(dummy);
    if( rank == 0 )
    {
        all_particles.clear();
        init_particles( n, particles );
        
        update_cells(particles, cells, n);
        for(int rankId = 0; rankId < n_proc; rankId++) 
        {
            partition_offsets[rankId] = all_particles.size();
            int first_row = min(  rankId     * rows_per_thread, num_cells);
            int last_row  = min( (rankId+1)  * rows_per_thread, num_cells);

            list<particle_t> tmp;
            tmp.clear();
            get_particles_from_rows(first_row, last_row, &tmp, cells);
            all_particles.insert(all_particles.end(), tmp.begin(), tmp.end());
           // int amount = tmp.size();
            partition_sizes[rankId] = tmp.size();

        }
       partition_offsets[n_proc] = n;
    }
    // broadcast all offsets ant sizes so we can scatter later.
    MPI_Bcast(partition_offsets, n_proc+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(partition_sizes, n_proc, MPI_INT, 0, MPI_COMM_WORLD);
    // get my_amount from the partition sizes array and rezise the my_particles vector.
    my_amount = partition_sizes[rank];
    my_particles.resize(my_amount);
    
    my_particles.push_back(dummy);
    MPI_Scatterv( &all_particles.front(), partition_sizes, partition_offsets, PARTICLE, &my_particles.front(), my_amount, PARTICLE, 0, MPI_COMM_WORLD );   
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        // 
        //  collect all global data locally (not good idea to do)
        //
        MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        update_cells(particles, cells, n);
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        for( int i = 0; i < nlocal; i++ )
        {
            local[i].ax = local[i].ay = 0;
            apply_force(&local[i], cells);

            // for (int j = 0; j < n; j++ )
            //     apply_force( local[i], particles[j] );
        }
        
        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}