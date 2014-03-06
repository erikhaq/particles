#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include <iostream>
#include <typeinfo>

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force_to_particle( particle_t &particle, particle_t &neighbor )
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
    {
        return;
    }
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
    
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

int get_num_cells()
{
    return (int) ceil(size/cutoff);
}

void init_cell_matrix(CellMatrix &cells)
{
    int num_cells = cells.size();
    
    for(int i = 0; i < num_cells ; i++)
    {
        cells[i].resize(num_cells);
    }
    
}

void update_cells(particle_t *particles, CellMatrix &cells, int n)
{
    int i, j;
    int num_cells = cells.size();
    
    for(i = 0; i < num_cells; i++)
    {
        for(j = 0; j < num_cells; j++)
        {
            cells[i][j].clear();
            cells[i][j].resize(0);
        }
    }


    for(i = 0; i < n; i++)
    {
        Point p = get_cell_index(particles[i]);
        cells[p.y][p.x].push_back(&particles[i]);
    }
    
}

void update_cells_only(int first_particle, int last_particle, particle_t *particles, CellMatrix &cells)
{
    for(int i = first_particle; i < last_particle; i++)
    {
        Point p = get_cell_index(particles[i]);
        cells[p.y][p.x].push_back(&particles[i]);
    }
}

void clear_cells(CellMatrix &cells)
{    
    int num_rows = cells.size();
    
    for(int i = 0; i < num_rows; i++)
    {
        int num_cols = cells[i].size();
        for(int j = 0; j < num_cols; j++)
        {
            cells[i][j].clear();
            cells[i][j].resize(0);

        }
    }
}

void clear_cells(int start_row, int end_row, CellMatrix& cells)
{ 
    for(int i = start_row; i < end_row; i++)
    {
        int num_cols = cells[i].size();
        for(int j = 0; j < num_cols; j++)
        {
            cells[i][j].clear();
            cells[i][j].resize(0);

        }
    }
}

Point get_cell_index(particle_t &particle)
{
    Point p = {floor(particle.x/cutoff), floor(particle.y/cutoff)};
    return p;
}

void apply_force(particle_t *particle, CellMatrix &cells)
{
    Point center = get_cell_index(*particle);
    int num_cells = get_num_cells();
    int r_start = center.y - 1;                         // row start
    r_start = clamp<int>(r_start, 0, center.y);         
    int r_end = center.y + 1;                           // row end
    r_end = clamp<int>(r_end, 0, num_cells-1);

    int c_start = center.x - 1;                         // column start
    c_start = clamp<int>(c_start, 0, center.x);
    int c_end = center.x + 1;                           // column end
    c_end = clamp<int>(c_end, 0, num_cells-1);

    for(int r = r_start; r <= r_end; r++)
    {
        for(int c = c_start; c <= c_end; c++)
        {       
            if(!cells[r][c].empty())
            {
                Particles neighbors = cells[r][c];
                int num_parts = neighbors.size();
                for(int i = 0; i < num_parts; i++)
                {       
                    particle_t *neighbor = cells[r][c][i];
                    if(particle == neighbor)
                    {
                        // we don't want to check ourselves
                        continue;
                    }
                    // apply the force of neighbor on current particle
                    apply_force_to_particle(*particle, *neighbor);
                }      
            }
        }
    }

}

void get_particles_from_rows(int start_row, int end_row, Particles *in, CellMatrix& cells)
 {
    for( int r = start_row; r < end_row; r++) 
        {
            int num_cols = cells[r].size();
            for(int c = 0; c < num_cols; c++)
            {
                Particles cell_particles = cells[r][c];
                int num_parts = cell_particles.size();
                for(int i = 0; i < num_parts; i++)
                {
                   // particle_t *curr_particle = cell_particles[i];
                    in->push_back(cell_particles[i]);
                }
            }
        }
 }

void print_point(Point p)
{
    printf("Point: (%d, %d)\n", p.y, p.x);
}

void print_particle(particle_t *particle)
{
    printf("Particle: (%.8f, %.8f)\n",particle->x, particle->y );
    
}

bool is_same(particle_t *p1, particle_t *p2)
{
    return (p1->x == p2->x && p1->y == p2->y);
}

