#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif

__kernel void matmult(__global double* A, __global double* B, __global double* C,__global size_t* d_nA)
{
  int row = get_global_id(0); //i
  int col = get_global_id(1); //j
  int n = *d_nA;
 
  // Valor resultado local calculado por el hilo
  double Cvalue = 0;
  
  // Si se pasa de tamaño acaba
  if( (row * col) < (n * n) ) {
    int index_a,index_b;
  
    // Multiplica las matrices, cada hilo calcula un elemento
    for (int k = 0; k < n; k++){      
      index_a = row+k*n;
      index_b = col*n+k;
      Cvalue += A[ index_a ] * B[index_b];
     }
     
     // Almacena el valor
     C[row+n*col] = Cvalue;
  }
}
