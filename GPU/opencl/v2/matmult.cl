#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif

// Thread block size
#define BLOCK_SIZE 16

__kernel void matmult(__global double* A, __global double* B, __global double* C,__global size_t* d_nA)
{
  int n = *d_nA;

  // Block index
  int bx = get_group_id(0);
  int by = get_group_id(1);

  // Thread index
  int tx = get_local_id(0);
  int ty = get_local_id(1);

  // Indice de la primera submatriz de A procesada por el bloque
  int aBegin = n * BLOCK_SIZE * by;
  // Indice de la ultima submatriz de A procesada por el bloque
  int aEnd   = aBegin + n - 1;

  int aStep  = BLOCK_SIZE;

  // Indice de la primera submatriz de B procesada por el bloque
  int bBegin = BLOCK_SIZE * bx;
  // Indice de la ultima submatriz de B procesada por el bloque
  int bStep  = BLOCK_SIZE * n;

  double Csub;

  for (int a = aBegin, b = bBegin;
    a <= aEnd;
    a += aStep, b += bStep) 
  {

    // Declaration de la meoria local de As para almacenar la submatriz A
    __local double As[BLOCK_SIZE][BLOCK_SIZE];

    // Declaration de la meoria local de Bs para almacenar la submatriz B
    __local double Bs[BLOCK_SIZE][BLOCK_SIZE];

    // Copiar matrices locales
    As[ty][tx] = A[a + n * ty + tx];
    Bs[ty][tx] = B[b + n * ty + tx];

    // Sincronizar para que nos aseguremos de que la matrices se han copiado
    barrier(CLK_LOCAL_MEM_FENCE);

    // Multiplicacion de las dos submatrices locales
    // Cada hilo multiplicará un elemento de la submatriz
    for (int k = 0; k < BLOCK_SIZE; ++k){
      Csub += As[ty][k] * Bs[k][tx];
    }

    // Sincronizar para asegurar que se ha acabado la operacion
    barrier(CLK_LOCAL_MEM_FENCE);

  }

  // Escribir el resultado en C, pero cada hilo escribe un resultado
  int c = n * BLOCK_SIZE * by + BLOCK_SIZE * bx;
  C[c + n * ty + tx] = Csub;

}
