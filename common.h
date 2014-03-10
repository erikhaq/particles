#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>
#include <list>
#include <math.h>
using namespace std;

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

//
// cell data structure
//
typedef vector<particle_t*> Particles;
typedef vector<particle_t> ParticleList;
typedef vector<vector<Particles> > CellMatrix;
typedef struct {
	int x, y;
} Point;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force_to_particle( particle_t &particle, particle_t &neighbor );
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

// erik and ludwig super functions
int get_num_cells();
void init_cell_matrix(CellMatrix&);
void update_cells(particle_t*, CellMatrix&, int);
void update_cells_only(int, int, particle_t*, CellMatrix&);
void add_particles_to_cells(ParticleList&, CellMatrix& );
Point get_cell_index(particle_t&);
void apply_force(particle_t*, CellMatrix&);
void print_point(Point);
void print_particle(particle_t*);
bool is_same(particle_t*, particle_t*);
void clear_cells(CellMatrix&);
void clear_cells(int, int, CellMatrix&);
void get_particles_from_rows(int, int, Particles*, CellMatrix&);
void get_particles_from_rows(int , int , ParticleList* , CellMatrix&);

template <typename T>
T clamp(T in, T min, T max)
{
	return std::min(std::max(in, min), max);
}



#endif
