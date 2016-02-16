#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

/* 
 * Método de Jacobi para la ecuación de Poisson
 *  
 *   Suponemos definida una malla de (N+1)x(M+1) puntos, donde los puntos
 *   de la frontera tienen definida una condición de contorno.
 *  
 *   Esta función resuelve el sistema Ax=b mediante el método iterativo
 *   estacionario de Jacobi. La matriz A no se almacena explícitamente y
 *   se aplica de forma implícita para cada punto de la malla. El vector
 *   x tiene dimensión NxM y representa la solución de la ecuación de Poisson
 *   en cada uno de los puntos de la malla. El vector b es también de NxM
 *   y contiene el término h^2*f.
 *
 *   Suponemos que las condiciones de contorno son igual a 0 en toda la
 *   frontera del dominio.
 */
void jacobi_poisson(int N,int M,double *x, double *b, int maxit,double tol, int me, int np)
{
  int i, j, k, conv;
  double *t, s;
  double s_global;
  int process;

  MPI_Status status;
  MPI_Request requestEnvio;

  t = (double*)calloc(N*M,sizeof(double));

  k = 0;
  conv = 0;

  while(!conv && k<maxit) {

    if(me>0){
       process = me-1;
       MPI_Isend(t+1*M+1,M-2,MPI_DOUBLE,process,0,MPI_COMM_WORLD,&requestEnvio);
       MPI_Recv(t+1,M-2,MPI_DOUBLE,process,0,MPI_COMM_WORLD,&status);
    }
    if (me<np-1) {
      process = me+1;
      MPI_Isend(t+(N-2)*M+1,M-2,MPI_DOUBLE,process,0,MPI_COMM_WORLD,&requestEnvio);
      MPI_Recv(t+(N-1)*M+1,M-2,MPI_DOUBLE,process,0,MPI_COMM_WORLD,&status);

    }
    /* interior de la malla */
    for(i=1; i<N-1; i++) {
      t[i*M+0] = (b[i*M+0] + x[(i+1)*M+0] + x[(i-1)*M+0] + x[i*M+(0+1)])/4.0;
      for(j=1; j<M-1; j++) {
        t[i*M+j] = (b[i*M+j] + x[(i+1)*M+j] + x[(i-1)*M+j] + x[i*M+(j+1)] + x[i*M+(j-1)])/4.0;
      }
      t[i*M+(M-1)] = (b[i*M+(M-1)] + x[(i+1)*M+(M-1)] + x[(i-1)*M+(M-1)] + x[i*M+((M-1)-1)])/4.0;
    }


    /* criterio de parada: ||x_{k}-x_{k+1}||<tol */
    s = 0.0;
    for(i=1; i<N-1; i++) {
      for(j=1; j<M-1; j++) {
        s += (x[i*M+j]-t[i*M+j])*(x[i*M+j]-t[i*M+j]);
      }
    }
    
    MPI_Allreduce(&s, &s_global, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

    conv = (sqrt(s_global)<tol);

    printf("Error en iteración %d: %g\n", k, sqrt(s_global));

    /* siguiente iteración */
    k = k+1;
    for(i=1; i<N-1; i++) {
      for(j=1; j<M-1; j++) {
        x[i*M+j] = t[i*M+j];
      }
    }

  }

  free(t);
}

int main(int argc, char **argv) 
{
  int i, j, N=50, M=50;
  double *x, *b, h=0.01, f=1.5;

  int mit=10000;
  double tol=0.000001;

  int np,me;
  int n;

  /*Iniciar MPI*/
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de N */
     if ((N = atoi(argv[1])) < 0) N = 50;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de M */
     if ((M = atoi(argv[2])) < 0) M = 1;
  }
  if (argc > 3) { /* El usuario ha indicado el valor de M */
     if ((mit = atoi(argv[3])) < 1) mit = 10000;
  }
  if (argc > 4) { /* El usuario ha indicado el valor de M */
     if ((tol = atoi(argv[4])) < 0) tol=1e-6;
  }

  for(i=0; i<np; i++) {
    n = N/np;
    if(i<N%np)n++;
    if(i==me) break;
  }


  n=n+2;
  M=M+2;

  /* Reserva de memoria */
  x = (double*)calloc(n*M,sizeof(double));
  b = (double*)calloc(n*M,sizeof(double));

  /* Inicializar datos */
  for(i=1; i<n-1; i++) {
    for(j=1; j<M-1; j++) {
      b[i*M+j] = h*h*f;  /* suponemos que la función f es constante en todo el dominio */
    }
  }

  /* Imprimir inicialización */
/*  for(i=0; i<n; i++) {
    for(j=0; j<M; j++) {
      printf("%g ", b[i*M+j]);
    }
    printf("\n");
  }
*/
  /* Resolución del sistema por el método de Jacobi */
  jacobi_poisson(n,M,x,b,mit,tol,me,np);
  
  /* Imprimir solución */
/*  MPI_Barrier(MPI_COMM_WORLD);
  for(i=1; i<n-1; i++) {
    printf("Proceso %d ", me);
    for(j=1; j<M-1; j++) {
      printf("%g ", x[i*M+j]);
    }
    printf("\n");
  }
*/
  free(x);
  free(b);
  MPI_Finalize();
  
  return 0;
}

