#include <stdio.h>
#include <stdlib.h>

/* 
 * Multiplicación de una matriz banda por un vector
 *  w = A*v, con A matriz cuadrada de dimensión N y ancho de banda b
 *  Algoritmo orientado a filas
 */
void matvec(int N,int b,double *A, double *v, double *w)
{
  int i, j, li, ls;

  for(i=0; i<N; i++) {
    w[i] = 0.0;
    li = i-b<0? 0: i-b;  /* limite inferior */
    ls = i+b>N-1? N-1: i+b;  /* limite superior */
    for(j=li; j<=ls; j++) {
      w[i] += A[i*N+j]*v[j];
    }
  }
}

int main(int argc, char **argv) 
{
  int i, j, N=50, b=1;
  double *A, *v, *w;

  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de N */
     if ((N = atoi(argv[1])) < 0) N = 50;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de b */
     if ((b = atoi(argv[2])) < 0) b = 1;
  }
  if (b>=N) { /* Valor de b incorrecto */
    printf("Error: ancho de banda excesivo, N=%d, b=%d\n", N, b);
    exit(1);
  }

  /* Reserva de memoria */
  A = (double*)calloc(N*N,sizeof(double));
  v = (double*)calloc(N,sizeof(double));
  w = (double*)calloc(N,sizeof(double));

  /* Inicializar datos */
  for(i=0; i<N; i++) A[i*N+i] = 2*b;
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      if (i!=j && abs(i-j)<=b) A[i*N+j] = -1.0;
    }
  }
  for(i=0; i<N; i++) v[i] = 1.0;

  struct timeval t1, t2;
  double elapsedTime;

  gettimeofday(&t1, NULL);
  
  /* Multiplicación de matrices */
  matvec(N,b,A,v,w);
  
  gettimeofday(&t2, NULL);
  // compute and print the elapsed time in millisec
  elapsedTime = (t2.tv_sec - t1.tv_sec)+((t2.tv_usec - t1.tv_usec)/(1000.0*1000.0));
  printf("TIME TAKEN : %f seconds\n",elapsedTime);

  /* Imprimir solución */
  for(i=0; i<N; i++) printf("w[%d] = %g\n", i, w[i]);

  free(A);
  free(v);
  free(w);

  return 0;
}

