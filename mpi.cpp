#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"

#define NUM_TOP 1
#define NUM_BOTTOM 2
#define BUFF_TOP 3
#define BUFF_BOTTOM 4
#define BUFF_SIZE 5
#define EXCHANGE 6

void send_overlapping_particles(CellMatrix &cells, int top_row_overlapp, int bottom_row_overlapp, int myRank, int n_proc, MPI_Comm comm, MPI_Datatype type);
void exchange_particles(ParticleList &particles, int myRank, int target, ParticleList &out, MPI_Datatype type ,MPI_Comm comm);
//
//  benchmarking program
//

void send_overlapping_particles(CellMatrix &cells, int top_row_overlapp, int bottom_row_overlapp, int myRank, int n_proc, MPI_Comm comm, MPI_Datatype type)
{
    int num_bottom, num_top;
    MPI_Status status;
    ParticleList myBuff, otherBuff;
    myBuff.clear();
    otherBuff.clear();
    if((myRank % 2) == 0) // all even processes sends their bottom row to the process below first, then top
    {
        if(myRank < n_proc-1) // we don't want the last process to send its bottom row
        {
            
            myBuff.clear();
            get_particles_from_rows(bottom_row_overlapp-1, bottom_row_overlapp, &myBuff, cells);
            exchange_particles(myBuff, myRank, myRank+1, otherBuff, type, comm);

            add_particles_to_cells(otherBuff, cells); // add the reference particles to the cells,
        }
        if(myRank > 0) // we don't want the root to send its top row anywhere.
        {
            myBuff.clear();
            get_particles_from_rows(top_row_overlapp, top_row_overlapp+1, &myBuff, cells);
            exchange_particles(myBuff, myRank, myRank-1, otherBuff, type, comm);
            add_particles_to_cells(otherBuff, cells); // add the reference particles to the cells,
        }
    }
    if(( myRank % 2 ) == 1)
    {
        myBuff.clear();
        get_particles_from_rows(top_row_overlapp, top_row_overlapp+1, &myBuff, cells);
        exchange_particles(myBuff, myRank, myRank-1, otherBuff, type, comm);
        
        add_particles_to_cells(otherBuff, cells); // add the reference particles to the cells,
        if(myRank < n_proc-1) // we don't want the last process to send its bottom row
        {
            myBuff.clear();
            get_particles_from_rows(bottom_row_overlapp-1, bottom_row_overlapp, &myBuff, cells);
            exchange_particles(myBuff, myRank, myRank+1, otherBuff, type, comm);

            add_particles_to_cells(otherBuff, cells); // add the reference particles to the cells,
        }

    }

}
void synchronize_new_particles(ParticleList &new_particles, ParticleList &out_of_bounds_top, ParticleList &out_of_bounds_bottom, int myRank, int n_proc, MPI_Comm comm, MPI_Datatype type)
{
    
    MPI_Status status;
    ParticleList tmp;
    
    /*
        All even processes start by communicating with the process below them
        and then they communicate with the process above them. This is to avoid 
        deadlock. 
    */
    if((myRank % 2) == 0) // all even processes sends their bottom row to the process below first, then top
    {
        if(myRank < n_proc-1) // we don't want the last process to send its bottom row
        {
            tmp.clear();
            exchange_particles(out_of_bounds_bottom, myRank, myRank+1, tmp, type, comm);
            new_particles.insert(new_particles.end(), tmp.begin(), tmp.end()); // append the new particles to out buffer

        }
        if(myRank > 0) // we don't want the root to send its top row anywhere.
        {
            tmp.clear();
            exchange_particles(out_of_bounds_top, myRank, myRank-1, tmp, type, comm);
            new_particles.insert(new_particles.end(), tmp.begin(), tmp.end()); // append the new particles to out buffer
        }
    }

    /*
        All odd processes start by communicating with the process above them
        and then they communicate with the process below them. This is the reverse
        of what the even processes do and is needed to avoid deadlock.
    */
    if(( myRank % 2 ) == 1)
    {
        tmp.clear();
        exchange_particles(out_of_bounds_top, myRank, myRank-1, tmp, type, comm);
        new_particles.insert(new_particles.end(), tmp.begin(), tmp.end()); // append the new particles to out buffer
        if(myRank < n_proc-1) // we don't want the last process to send its bottom row
        {
            tmp.clear();
            exchange_particles(out_of_bounds_bottom, myRank, myRank+1, tmp, type, comm);
            new_particles.insert(new_particles.end(), tmp.begin(), tmp.end()); // append the new particles to out buffer
        }

    }
}
void exchange_particles(ParticleList &particles, int myRank, int target, ParticleList &out, MPI_Datatype type ,MPI_Comm comm)
{
    MPI_Status status;
    int otherSize, mySize;
    mySize = particles.size();
    out.clear();
    MPI_Sendrecv(&mySize, 1, MPI_INT, target, BUFF_SIZE, &otherSize, 1, MPI_INT, target, BUFF_SIZE, comm, &status);
    out.resize(otherSize);
    // printf("%d sent %d and got %d from %d \n", myRank, mySize, otherSize, target);
    MPI_Sendrecv(&particles.front(), mySize, type, target, EXCHANGE, &out.front(), otherSize, type, target, EXCHANGE,comm, &status );

}
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
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE ); // define the datatype PARTICLE as struct of 6 doubles
    MPI_Type_commit( &PARTICLE );
    


    CellMatrix cells(num_cells);
   

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
    ParticleList my_particles;
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)

    ParticleList all_particles;
    all_particles.clear();
    init_cell_matrix(cells);
   
    if( rank == 0 )
    {
        init_particles( n, particles );
        
        update_cells(particles, cells, n);
        for(int rankId = 0; rankId < n_proc; rankId++) 
        {
            partition_offsets[rankId] = all_particles.size();
            int first_row = min(  rankId     * rows_per_thread, num_cells);
            int last_row  = min( (rankId+1)  * rows_per_thread, num_cells);

            ParticleList tmp;
            tmp.clear();
            get_particles_from_rows(first_row, last_row, &tmp, cells);
            all_particles.insert(all_particles.end(), tmp.begin(), tmp.end());
            int amount = tmp.size();
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
  
    MPI_Scatterv( &all_particles.front(), partition_sizes, partition_offsets, PARTICLE, &my_particles.front(), my_amount, PARTICLE, 0, MPI_COMM_WORLD );    
   /* clears all the rows that the process is responsible for and also the overlapping
      rows that this process needs for its computatuins. The clear cells function clamps 
      the row values so it does not go out of bounds since the first and last process does 
      not have a overlapping top and bottom row respectivly.
    */
    clear_cells(top_row-1, bottom_row+1, cells); 
    add_particles_to_cells(my_particles, cells);

    
    send_overlapping_particles(cells, top_row, bottom_row, rank, n_proc, MPI_COMM_WORLD, PARTICLE);

    // //
    // //  simulate a number of time steps
    // //
      double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        // 
        //  collect all global data locally (not good idea to do)
        //
        
        MPI_Barrier(MPI_COMM_WORLD); // wait in order to synchronize with same frame.


       // MPI_Gather(&my_particles.front(), my_particles.size(), PARTICLE, particles, n, PARTICLE, 0, MPI_COMM_WORLD);


        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        ParticleList::iterator iter = my_particles.begin();
        while(iter != my_particles.end())
        {
            particle_t *curr_particle = &(*iter);
            curr_particle->ax = 0;
            curr_particle->ay = 0;
            // local[i].ax = local[i].ay = 0;
            //  apply_force(&local[i], cells);
             apply_force(curr_particle, cells);
             ++iter;

        }
        
        //
        //  move particles
        
        ParticleList out_of_bounds_top;
        ParticleList out_of_bounds_bottom;
        out_of_bounds_top.clear();
        out_of_bounds_bottom.clear();
        iter = my_particles.begin();
        // printf("%d will start moving particles\n", rank);
        while(iter != my_particles.end())
        {
            
           // particle_t *curr_particle = &(*iter);
            move(*iter);
            Point p = get_cell_index(*iter);
            if(p.y < top_row ) // the particle moved out of our bounds
            {
                
                particle_t tmp = (*iter);
                out_of_bounds_top.push_back(tmp);
                // my_particles.erase(iter++);
                iter = my_particles.erase(iter);
                // printf("bounds top after %d\n", out_of_bounds_top.size());
            }
            else if(p.y >= bottom_row)
            {
                // printf("%d :outofboundds bottom\n", rank);
                out_of_bounds_bottom.push_back(*iter);
                // my_particles.erase(iter++);
                iter = my_particles.erase(iter);
            }
            else
            {
                ++iter;
            }
            // ++iter;

        }
        // printf("%d is done moving particles\n", rank);
        // clear_cells(top_row-1, bottom_row+1, cells); 
        // ParticleList new_particles;
        // new_particles.clear();
        // // printf("Hejsan %d\n", rank);
        // synchronize_new_particles(new_particles, out_of_bounds_top, out_of_bounds_bottom, rank, n_proc, MPI_COMM_WORLD, PARTICLE);
        //  // printf("Hejsan2 %d\n", rank);
        // MPI_Barrier(MPI_COMM_WORLD);
        // my_particles.insert(my_particles.end(), new_particles.begin(), new_particles.end());
        // add_particles_to_cells(my_particles, cells);
        // send_overlapping_particles(cells, top_row, bottom_row, rank, n_proc, MPI_COMM_WORLD, PARTICLE);
        // for( int i = 0; i < my_amount; i++ )
        // {
        //     ParticleList out_of_bounds;
        //     out_of_bounds.clear();
        //     particle_t curr = my_particles[i];
        //     move( curr );
        //     Point p = get_cell_index(curr);
        //     if(p.y < first_row || p.y >= last_row) // the particle moved out of our bounds

        // }

        
    }
    simulation_time = read_timer( ) - simulation_time;
     MPI_Barrier(MPI_COMM_WORLD);
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    free( partition_offsets );
    free( partition_sizes );
    // free( local );

    // if(rank == 0)
    free( particles );
    

    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
