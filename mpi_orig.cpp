#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <list>
#include "common.h"

using namespace std;
#define NEW 1
#define BOUND 2
#define DONE 3
#define SAVE 4
//
//  benchmarking program
//

// void add_particles_mpi(list<particle_t> &particles, CellMatrix &cells)
// {
//     list<particle_t>::iterator iter = particles.begin();
//     while(iter != particles.end())
//     {
//         particle_t tmp = (*iter);
//         Point p = get_cell_index(tmp);
//         cells[p.y][p.x].push_back(&(*iter));
//         ++iter;
//     }
// }
void add_particles_mpi(ParticleList &particles, CellMatrix &cells)
{
    ParticleList::iterator iter = particles.begin();
    while(iter != particles.end())
    {
        particle_t tmp = (*iter);
        Point p = get_cell_index(tmp);
        cells[p.y][p.x].push_back(&(*iter));
        ++iter;
    }
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
    // FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
     FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    
    // list<particle_t> reference_particles;
    // list<particle_t> my_particles;
    ParticleList reference_particles;
    ParticleList my_particles;

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

    ParticleList all_particles;
    all_particles.clear();
    init_cell_matrix(cells);
    // all_particles.push_back(dummy);
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

            ParticleList tmp;
            tmp.clear();
            get_particles_from_rows(first_row, last_row, &tmp, cells);
           // int amount = tmp.size();
            all_particles.insert(all_particles.end(), tmp.begin(), tmp.end());
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
    
    // my_particles.push_back(dummy);
    MPI_Scatterv( &all_particles.front(), partition_sizes, partition_offsets, PARTICLE, &my_particles.front(), my_amount, PARTICLE, 0, MPI_COMM_WORLD );   
    // printf("r: %d - size: %d\n",rank,  my_particles.size());
    
    clear_cells(top_row-1, bottom_row+1, cells); 
    add_particles_mpi(my_particles, cells);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        // 
        //  collect all global data locally (not good idea to do)
        //

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( fsave && (step%SAVEFREQ) == 0 )
        {


            if(rank == 0)
            {
                int num_waiting = n - my_particles.size();
                int i;
                int count = num_waiting;
                // MPI_Request *requests =(MPI_Request*) malloc(num_waiting * sizeof(MPI_Request));
                for(i = 0; i < my_particles.size(); i++)
                {
                    particles[i] = my_particles[i];
                }
                particle_t other;
                while(num_waiting > 0)
                {
                    MPI_Status stat;
                    MPI_Recv(&other, 1, PARTICLE, MPI_ANY_SOURCE, SAVE, MPI_COMM_WORLD, &stat);
                    particles[i] = other;
                    // MPI_Irecv(&particles[i], 1, PARTICLE, MPI_ANY_SOURCE, SAVE, MPI_COMM_WORLD, &requests[num_waiting-1]);
                    i++;
                    num_waiting--;
                    // printf("num_waiting: %d\n", num_waiting);
                }
                // MPI_Waitall(count,requests, MPI_STATUSES_IGNORE);
                save( fsave, n, particles );
            }
            else 
            {
                for(int i = 0; i  < my_particles.size(); i++)
                {
                    MPI_Request req;
                    MPI_Isend(&my_particles[i], 1 , PARTICLE, 0, SAVE, MPI_COMM_WORLD, &req );
                }
            }

            
        }
        
        //
        //  compute all forces
        //
        // list<particle_t>::iterator iter = my_particles.begin();
        ParticleList::iterator iter = my_particles.begin();
        while(iter != my_particles.end())
        {
            
            particle_t *curr_particle = &(*iter);
            curr_particle->ax = 0;
            curr_particle->ay = 0;
            
             apply_force(curr_particle, cells);
             ++iter;

        }
        reference_particles.clear();
        ParticleList out_of_bounds;
        out_of_bounds.clear();
        ParticleList bounds;
        bounds.clear();
        iter = my_particles.begin();
        while(iter != my_particles.end())
        {
            
           // particle_t *curr_particle = &(*iter);
            move(*iter);
            Point p = get_cell_index(*iter);
            if(p.y < top_row || p.y >= bottom_row)
            {
                particle_t tmp = (*iter);
                int index = out_of_bounds.size();
                out_of_bounds.push_back(tmp);
                iter = my_particles.erase(iter);

                int target = (p.y < top_row)? rank-1 : rank+1;
                MPI_Request request;
                MPI_Isend(&out_of_bounds[index], 1, PARTICLE, target, NEW, MPI_COMM_WORLD, &request ); // non blocking send
                continue;
            }
            else if(p.y == top_row && top_row > 0) // send our top row particles to the process above except if we are root
            {
                
            
                particle_t tmp = (*iter);
                int index = bounds.size();
                bounds.push_back(tmp);

                
                MPI_Request request;
                MPI_Isend(&bounds[index], 1, PARTICLE, rank-1, BOUND, MPI_COMM_WORLD, &request ); // non blocking send
            }
            else if(p.y == bottom_row-1 && rank < n_proc-1) // send our
            {
                particle_t tmp = (*iter);
                int index = bounds.size();
                bounds.push_back(tmp);

                
                MPI_Request request;
                MPI_Isend(&bounds[index], 1, PARTICLE, rank+1, BOUND, MPI_COMM_WORLD, &request ); // non blocking send
            }

            ++iter;
             
        }
        MPI_Request req;
        //
        if(top_row > 0) MPI_Isend(&dummy, 1, PARTICLE, rank - 1, DONE, MPI_COMM_WORLD, &req);
        if(bottom_row < num_cells) MPI_Isend(&dummy, 1, PARTICLE, rank + 1, DONE, MPI_COMM_WORLD, &req);

        int isDone = 2;
        if(top_row  == 0) isDone--;
        if(bottom_row == num_cells) isDone--;

        particle_t new_particle;
        while(isDone > 0)
        {
            MPI_Status status;
            MPI_Recv(&new_particle, 1, PARTICLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if(status.MPI_TAG == DONE)
            {
                isDone--;
                continue;
            }
            else if(status.MPI_TAG == NEW)
            {
                my_particles.push_back(new_particle);
            }
            else if(status.MPI_TAG == BOUND)
            {
                reference_particles.push_back(new_particle);
            }

        }
    MPI_Barrier(MPI_COMM_WORLD); // wait in order to synchronize with same frame.   
    clear_cells(top_row-1, bottom_row+1, cells); 
    add_particles_mpi(my_particles, cells);
    add_particles_mpi(reference_particles, cells);

  
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    free( partition_offsets );
    free( partition_sizes );
    //free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}