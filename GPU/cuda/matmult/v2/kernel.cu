__global__ void matmultKernel(double* A, double* B, double* C, int n) {
	
	// Bloque al que pertenecen los threads
	int bx = blockIdx.x;
	int by = blockIdx.y;
	
	// Identificador del Thread dentro del bloque
	int tx = threadIdx.x; 
	int ty = threadIdx.y;

	int row = by*blockDim.y+ty;
	int col = bx*blockDim.x+tx;

	// Valor resultado local calculado por el hilo
	double Cvalue = 0;
	
	// Si se pasa de tamaño acaba
	//if(row < n && col < n) {
	if(row >= n || col >= n) return;

		int index_a,index_b;
	
		// Multiplica las matrices, cada hilo calcula un elemento
		for (int k = 0; k < n; k++){			
			index_a = col+k*n;
			index_b = row*n+k;
			Cvalue += A[ index_a ] * B[index_b];
			//printf("A: %f\t B: %f \t C: %f\n",A[index_a],B[index_b],Cvalue); 
	   }
	   
	   // Almacena el valor
	   C[row*n+col] = Cvalue;
	//}
}