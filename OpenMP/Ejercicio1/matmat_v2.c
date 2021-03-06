#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/* 
 * Multiplicación de matrices
 *  C = A*B, con A, B, C matrices cuadradas de dimensión n
 */
void matmat(int n,double *A, double *B, double *C)
{
  int i, j, k;

  for(i=0; i<n; i++) {
    #pragma omp parallel for private(k)
    for(j=0; j<n; j++) {
      C[i*n+j] = 0.0;
      for(k=0; k<n; k++) {
        C[i*n+j] += A[i*n+k]*B[k*n+j];
      }
    }
  }
}

int main(int argc, char **argv) 
{
  int i, j, n=50;
  double *A, *B, *C, time;

  /* Extracción de argumentos */
  if (argc == 2) { /* El usuario ha indicado el valor de n */
     if ((n = atoi(argv[1])) < 0) n = 50;
  }

  /* Creación de las matrices */
  A = (double*)malloc(n*n*sizeof(double));
  B = (double*)malloc(n*n*sizeof(double));
  C = (double*)malloc(n*n*sizeof(double));

  /* Inicializar matrices */
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      A[i+n*j] = drand48();
      B[i+n*j] = drand48();
    }
  }

  /* Multiplicación de matrices */
  time = omp_get_wtime();
  matmat(n,A,B,C);
  time = omp_get_wtime() - time;  
  printf("%d;%d;%f\n", omp_get_max_threads(), n, time);

  free(A);
  free(B);
  free(C);

  return 0;
}

