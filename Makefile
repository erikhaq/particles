#
# Computers with Red Hat Enterprise Linux 5 in the computer room 648, KTH Forum, Kista
#

CC = g++-4.9
# MPCC =  mpicc -cc=g++44
# MPCC =  mpicc -lmpi -lmpi_cxx
MPCC =  mpicxx
OPENMP = -fopenmp
LIBS = -lm
CFLAGS = -O3

TARGETS = serial pthreads openmp mpi mpi2
# TARGETS = serial pthreads openmp mpi mpi_orig

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o
mpi_orig: mpi_orig.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi_orig.o common.o
mpi2: mpi2.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi2.o common.o

openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
mpi_orig.o: mpi_orig.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi_orig.cpp	
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.out
