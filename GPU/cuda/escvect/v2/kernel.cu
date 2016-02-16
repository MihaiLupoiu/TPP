/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

__global__ void escvectKernel(float* A, float* B, float* result, int n) {

	float res = 0;
	int i = blockDim.x*blockIdx.x+threadIdx.x;
	__shared__ float sdata[1];


	if (threadIdx.x==0) {
    	sdata[0]=0;
	}
	 
	if(i<n){
		res = A[i]*B[i];
		//printf("A[%d] = %f, B[%d] = %f, res = %f\n",i,A[i],i,B[i],res);
		__syncthreads();
		atomicAdd(&sdata[0],res);
		__syncthreads();
	} 

	
	if (threadIdx.x==0) {
		atomicAdd(result,sdata[0]);
	}
	
}