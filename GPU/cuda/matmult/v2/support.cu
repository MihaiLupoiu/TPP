/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "support.h"

void verify(double *A, double *B, double *C, int n) {

  const double relativeTolerance = 1e-6;
  
  double sum;

  for (int i = 0; i < n; ++i) {
  	for (int j = 0; j < n; j++) {
		sum = 0;
        	for (int k = 0; k < n; k++) {
                	sum = sum + A[i+k*n]*B[k+j*n];
            	}
		
			double relativeError = (sum - C[i+j*n])/sum;
	            	
			if (relativeError > relativeTolerance || relativeError < -relativeTolerance) {
	                	printf("TEST FAILED\n\n");
	                	printf("Error %f \t C[%d][%d] = %f \t sum= %f \n\n",relativeError,i,j,C[i*n+j],sum);
	                    exit(0);
			}
		}
  }
  printf("TEST PASSED\n\n");
  
}

void startTime(Timer* timer) {
    gettimeofday(&(timer->startTime), NULL);
}

void stopTime(Timer* timer) {
    gettimeofday(&(timer->endTime), NULL);
}

float elapsedTime(Timer timer) {
    return ((float) ((timer.endTime.tv_sec - timer.startTime.tv_sec) \
                + (timer.endTime.tv_usec - timer.startTime.tv_usec)/1.0e6));
}

