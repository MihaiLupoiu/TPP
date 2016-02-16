/**
 * transpose of 1-d array.
 * E,g: (1 2 3 4 5 6 7 8 9 10) => (10 9 8 7 6 5 4 3 2 1)
 */
#define   NDIM         1
#define   TOTALELEMS   197
#define   MAXDIM       7
#define   MAXPROC      128

#include <stdio.h>
#include <math.h>
#include "ga.h"
#include "macdecls.h"

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif

void verify(int g_a, int g_b);

void TRANSPOSE1D() {
    
    int dims[MAXDIM], chunk[MAXDIM], ld[MAXDIM], lo[MAXDIM], hi[MAXDIM];
    int lo1[MAXDIM], hi1[MAXDIM], lo2[MAXDIM], hi2[MAXDIM];
    int g_a, g_b, a[MAXPROC*TOTALELEMS],b[MAXPROC*TOTALELEMS];
    int nelem, i;    
    int me, nprocs;

    ga_nbhdl_t get_handle, put_handle;

    /* Find local processor ID and number of processors */
    /* ### assign the local processor ID to the int variable "me"
     * ### and the total number of processors to the int variable
     * ### "nprocs" */
    
    me     = GA_Nodeid();
    nprocs = GA_Nnodes();

    /* Configure array dimensions. Force an unequal data distribution */
    dims[0]  = nprocs*TOTALELEMS + nprocs/2;
    ld[0]    = dims[0];
    chunk[0] = TOTALELEMS; /* minimum data on each process */
 
    /* create a global array g_a and duplicate it to get g_b */
    /* ### create GA of integers with dimension "NDIM" and size "dims" with
     * ### minimum block size "chunk" and assign the handle to the
     * ### integer variable "g_a". Then create a second global array
     * ### assigned to the integer handle "g_b" by duplicating "g_a".
     * ### Assign the names "Array A" and "Array B" to "g_a" and "g_b". */

    g_a=NGA_Create(C_INT, 1, dims, "array A", chunk);
    g_b=GA_Duplicate(g_a,"array B");

    if (!g_a) GA_Error("create failed: A", NDIM);
    if (me==0) printf("  Created Array A\n");
    
    if (! g_b) GA_Error("duplicate failed",NDIM);
    if (me==0) printf("  Created Array B\n");
 
    /* initialize data in g_a */
    if (me==0) {
       printf("  Initializing matrix A\n");
       for(i=0; i<dims[0]; i++) a[i] = i;
       lo[0]  = 0;
       hi[0] = dims[0]-1;
     /* ### copy the contents of array "a" into the portion of global array
      * ### "g_a" described by "lo" and "hi". Use the array of strides
      * ### "ld" to describe the physical layout of array "a". */

      NGA_NbPut(g_a,lo,hi,a,ld,&put_handle);

    }

    /* Synchronize all processors to guarantee that everyone has data
       before proceeding to the next step. */

    /* ### synchronize all processors */
    GA_Sync();

    /* Start initial phase of inversion by inverting the data held locally on
       each processor. Start by finding out which data each processor owns. */

    /* ### find out which block of data my node ("me") owns for the global
     * ### array "g_a" and store the contents in the integer arrays "lo1" and
     * ### "hi1". */

    NGA_Distribution(g_a,me,lo1,hi1);

    /* Get locally held data and copy it into local buffer a  */

    /* ### use the arrays "lo1" and "hi1" to copy the locally held block of data
     * ### from the global array "g_a" into the local array "a". Use the array
     * ### of strides "ld" to describe the physical layout of "a". */

    NGA_NbGet(g_a,lo1,hi1,a,ld,&get_handle);
    NGA_NbWait(&get_handle);


    /* Invert data locally */
    nelem = hi1[0] - lo1[0] + 1;
    for (i=0; i<nelem; i++) b[i] = a[nelem-1-i];
    
    /* Invert data globally by copying locally inverted blocks into
     * their inverted positions in the GA */
    lo2[0] = dims[0] - hi1[0] -1;
    hi2[0] = dims[0] - lo1[0] -1;

    /* ### copy data from the local array "b" into the block of the global
     * ### array "g_b" described by the integer arrays "lo2" and "hi2". Use
     * ### the array of strides "ld" to describe the physical layout of "b". */

    NGA_NbPut(g_b,lo2,hi2,b,ld,&put_handle);

    /* Synchronize all processors to make sure inversion is complete */
    /* ### synchronize all processors */
    
    GA_Sync();

    /* Check to see if inversion is correct */
    if(me == 0) verify(g_a, g_b);
    
    /* Deallocate arrays */
    /* ### destroy global arrays "g_a" and "g_b" */
    GA_Destroy(g_a);
    GA_Destroy(g_b);
}

/*
 * Check to see if inversion is correct. Start by copying g_a into local
 * buffer a, and g_b into local buffer b.
 */
void verify(int g_a, int g_b) {

    int i, type, ndim, dims[MAXDIM], lo[MAXDIM], hi[MAXDIM], ld[MAXDIM];
    int a[MAXPROC*TOTALELEMS],b[MAXPROC*TOTALELEMS];
    ga_nbhdl_t get_handle;
    
    /* Get dimensions of GA */
    NGA_Inquire(g_a, &type, &ndim, dims);

    lo[0] = 0;
    hi[0] = dims[0]-1;
    /* ### copy the block of data described by the arrays "lo" and "hi" from
     * ### the global array "g_a" into the local array "a". Copy the same block
     * ### of data from "g_b" into the local array "b". Use the array of strides
     * ### "ld" to describe the physical layout of "a" and "b". */

    NGA_NbGet(g_a,lo,hi,a,ld,&get_handle);
    NGA_NbWait(&get_handle);
    NGA_NbGet(g_b,lo,hi,b,ld,&get_handle);
    NGA_NbWait(&get_handle);

    for(i=0; i<dims[0]; i++)
       if (a[i] != b[dims[0]-i-1]) 
       {
          printf("Mismatch: a[%d]=%d is not equal to b[%d]=%d\n",
                 i, a[i], dims[0]-i-1, b[dims[0]-i-1]);
          GA_Error("verify failed",0);
       }
    
    printf("  Transpose OK\n");
}


int main(int argc, char **argv) {
    int heap=300000, stack=300000;
    int me, nprocs;
    
    /* Step1: Initialize Message Passing library */
#ifdef MPI
    MPI_Init(&argc, &argv);   /* initialize MPI */
#else
    PBEGIN_(argc, argv);      /* initialize TCGMSG */
#endif

    /* Step2: Initialize GA */
    /* ### intialize the GA library */
    GA_Initialize();

    
    /* Step3: Initialize Memory Allocator (MA) */
    if(! MA_init(C_DBL, stack, heap) ) GA_Error("MA_init failed",stack+heap);

    /* ### assign the local processor ID to the int variable "me"
     * ### and the total number of processors to the int variable
     * ### "nprocs" */
    me     = GA_Nodeid();
    nprocs = GA_Nnodes();
    if(me==0) {
       printf("\nUsing %d processes\n\n", nprocs); fflush(stdout);
    }
    
       
    TRANSPOSE1D();
    
    if(me==0)printf("\nTerminating ..\n");
    /* ### terminate the GA library */
    GA_Terminate();

#ifdef MPI
    MPI_Finalize();    
#else
    PEND_();
#endif
}
