/************************************************************
This is a simple hello world program. Each processor prints out 
it's rank and the size of the current MPI run (Total number of
processors).
************************************************************/
#include <stdio.h>
#include "mpi.h"
 
main(int argc, char* argv[]) {
  int myid;       /* rank of the proces */
  int numprocs;   /* nb of process */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  /* print out my rank and this run's PE size */
  printf("I am process %d from %d \n",myid,numprocs);
  MPI_Finalize();  
}

