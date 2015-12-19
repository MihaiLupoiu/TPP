/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#include <stdio.h>
#include "support.h"
#include "kernel.cu"

void printMatrixCustom(float* matrix,size_t size_r,size_t size_c){
    printf("\n *************** MATRIX ****************\n\n");
    int i,j;
    for(i = 0; i < size_r; i++) {
        for (j = 0; j < size_c; ++j) {
            printf(" %f ",matrix[i*size_c+j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char**argv) {

    Timer timer;
    cudaError_t cuda_ret;

    // Initialize host variables ----------------------------------------------

    printf("\nSetting up the problem..."); fflush(stdout);
    startTime(&timer);

    unsigned int n;
    if(argc == 1) {
        n = 1000;
    } else if(argc == 2) {
        n = atoi(argv[1]);
    } else {
        printf("\n    Invalid input parameters!"
           "\n    Usage: ./matmult               # Matrix of size 1,000 x 1,000 is used"
           "\n    Usage: ./matmult <m>           # Matrix of size m x m is used"
           "\n");
        exit(0);
    }

    float* A_h = (float*) malloc( sizeof(float)*n*n );
    for (unsigned int i=0; i < n*n; i++) { A_h[i] = (rand()%100)/100.00; }

    float* B_h = (float*) malloc( sizeof(float)*n*n );
    for (unsigned int i=0; i < n*n; i++) { B_h[i] = (rand()%100)/100.00; }

    float* C_h = (float*) malloc( sizeof(float)*n*n );

    stopTime(&timer); printf("%f s\n", elapsedTime(timer));
    printf("    Matrix size = %u x %u\n", n, n);
    
    printMatrixCustom(A_h,n,n);
    printMatrixCustom(B_h,n,n);

    // Allocate device variables ----------------------------------------------

    printf("Allocating device variables..."); fflush(stdout);
    startTime(&timer);

    //INSERT CODE HERE

	int size = n*n*sizeof(float);

	float* A_d;
	cudaMalloc((void**)&A_d,size);

	float* B_d;
	cudaMalloc((void**)&B_d,size);

	float* C_d;
	cudaMalloc((void**)&C_d,size);

    cudaDeviceSynchronize();
    stopTime(&timer); printf("%f s\n", elapsedTime(timer));

    // Copy host variables to device ------------------------------------------

    printf("Copying data from host to device..."); fflush(stdout);
    startTime(&timer);

    //INSERT CODE HERE

	cudaMemcpy(A_d, A_h, size, cudaMemcpyHostToDevice);
	cudaMemcpy(B_d, B_h, size, cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    stopTime(&timer); printf("%f s\n", elapsedTime(timer));

    // Launch kernel ----------------------------------------------------------

    printf("Launching kernel..."); fflush(stdout);
    startTime(&timer);

    //INSERT CODE HERE

	int nhilos = 4;
	int nbloques = ceil(float(n)/nhilos);

	int BLOCK_SIZE = nhilos;

	// llamada al kernel
	dim3 thrds(BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid(n/thrds.x, n/thrds.y);
	
	matmultKernel<<<grid,thrds>>>(A_d, B_d, C_d,n,n);


    cuda_ret = cudaDeviceSynchronize();
    if(cuda_ret != cudaSuccess) FATAL("Unable to launch kernel");
    stopTime(&timer); printf("%f s\n", elapsedTime(timer));

    // Copy device variables from host ----------------------------------------

    printf("Copying data from device to host..."); fflush(stdout);
    startTime(&timer);

    //INSERT CODE HERE

	cudaMemcpy(C_h, C_d, size, cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    stopTime(&timer); printf("%f s\n", elapsedTime(timer));

    // Verify correctness -----------------------------------------------------

    printMatrixCustom(C_h,n,n);

    printf("Verifying results..."); fflush(stdout);

    verify(A_h, B_h, C_h, n);

    // Free memory ------------------------------------------------------------

    free(A_h);
    free(B_h);
    free(C_h);

    //INSERT CODE HERE

	cudaFree(A_d);
	cudaFree(B_d);
	cudaFree(C_d);


    return 0;

}

