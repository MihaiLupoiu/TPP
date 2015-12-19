/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

__global__ void matmultKernel(float* A, float* B, float* C, int nA, int nB) {

	//int i = blockDim.x*blockIdx.x+threadIdx.x;
	//if(i<n)	C[i] = A[i]+B[i];
	
	//Bloque al que pertenecen los threads
	int bx = blockIdx.x;
	int by = blockIdx.y;
	
	//Identificador del Thread dentro del bloque
	int tx = threadIdx.x; 
	int ty = threadIdx.y;
	
	int BLOCK_SIZE = blockDim.x;
	
	int aBegin = nA*BLOCK_SIZE*by;     // Indice de primera submatriz
	int bBegin = BLOCK_SIZE*bx;

	// Valor resultado local calculado por el hilo
	float Csub = 0;
	
	int index_a, index_b;
	
	//printf("Hola!\n");
	
	// Recorrer todas las submatrices de A y B
	for (int a = aBegin, b = bBegin; a <= aBegin + nA - 1; a += BLOCK_SIZE, b += BLOCK_SIZE*nB) {			
		
		// Multiplica las matrices, cada hilo calcula un elemento
		for (int k = 0; k < BLOCK_SIZE; k++){			
			
			index_a = (a + nA * ty + tx)+k;
			index_b = (b + nB * tx + ty)+k;
			
        	Csub += A[ index_a ] * B[index_b];
        	
			printf("A: %f\t B: %f \t C: %f\n",A[index_a],B[index_b],Csub);
        	
        }
   }
   
   // Almacena sub-matriz calculada, cada hilo copia un elemento
   int c = nB * BLOCK_SIZE * by + BLOCK_SIZE * bx;
   C[c + nB * ty + tx] = Csub;


}

