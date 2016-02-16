#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/* Variables globales:
 *   np: número de procesos
 *   me: identificador de proceso
 */
int np, me;

/* 
 * Multiplicación PARALELA de una matriz banda por un vector
 *  w = A*v, con A matriz cuadrada de dimensión N y ancho de banda b
 *  Algoritmo orientado a filas
 *
 *  La matriz y los vectores estan distribuidos de la siguiente forma:
 *  - Cada procesador tiene n elementos consecutivos de los vectores,
 *    empezando en la posicion p
 *  - Cada procesador tiene n filas completas de la matriz, desde la
 *    fila p
 *  - A representa la porcion local de la matriz
 *  - v es un vector de longitud n+2*b con la siguiente estuctura:
 *      | b elem recib de P_{i-1} | n elem locales | b elem recib de P_{i+1} |
 *  - w es un vector con la misma estructura que v
 */
void parmatvec(int N,int n,int p,int b,double *A, double *v, double *w)
{
  int i, j, li, ls;
  MPI_Status status;
  MPI_Request requestEnvio;

  /* comunicación con vecino superior */
  if (me>0) {
    MPI_Isend(v,b,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,&requestEnvio);
  }
  
  for(i=0; i<n; i++) {
    w[i] = 0.0;
    li = p+i;  /* limite inferior */ // para que siempre empiece en p o en el lugar que le toque.
    ls = i+p+b>p+n-1? p+n-1: i+p+b;  /* limite superior */ //para acabar en p+n-1
    for(j=li; j<=ls; j++) {
      w[i] += A[i*N+j]*v[j-p];
    }
  }

  for(i=1; i<n; i++) {
    li = i+p-b<p? p: i+p-b;  /* limite inferior */
    ls = p+i-1;  /* limite superior */
    for(j=li; j<=ls; j++) {
      w[i] += A[(j-p)*N+(i+p)]*v[j-p];
    }
  }  
  
  /* comunicación con vecino inferior */
  if (me<np-1) {
    MPI_Recv(v+n,b,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,&status);
//    MPI_Wait(&requestEnvio,&status); Irecv
  }


  if (me != np-1) {   
	  
	  for(i=n-b; i<n; i++) {
	    li = p+n;  /* limite inferior */
	    ls = i+p+b;  /* limite superior */
	    for(j=li; j<=ls; j++) {
	      w[i] += A[i*N+j]*v[j-p];
	    }
	  }
	
	  for(i=n; i<n+b; i++) {
	    li = i+p-b;  /* limite inferior */
	    ls = p+n-1;  /* limite superior */
	    for(j=li; j<=ls; j++) {
	      w[i] += A[(j-p)*N+(i+p)]*v[j-p];
	    }
	  }
	  
	  /* comunicación con vecino superior */
	  MPI_Isend(w+n,b,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,&requestEnvio);
  }
  
  /* comunicación con vecino inferior */
  if (me>0) {
    MPI_Recv(v,b,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,&status);
//    MPI_Wait(&requestEnvio,&status); Irecv

    for(i=0; i<b; i++) {
    	w[i]+=v[i];
    }
  }  
  
}

int main(int argc, char **argv) 
{
  int i, j, n, N=50, b=1, p;
  double *A, *v, *w;

  /* Iniciar MPI */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
	
  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de n */
     if ((N = atoi(argv[1])) < 0) N = 50;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de b */
     if ((b = atoi(argv[2])) < 0) b = 1;
  }
  if (b>=N/np) { /* Valor de b incorrecto */
    printf("Error: ancho de banda excesivo, N/np=%d, b=%d\n", N/np, b);
    MPI_Finalize();
    exit(1);
  }

  p = 0;
  for (i=0; i<np; i++) {
    n = N/np;
    if (i<N%np) n++;
    if (i==me) break;
    p += n;
  }
  printf("[Proc %d] tamaño local: n=%d, fila inicial: p=%d\n", me, n, p);

  /* Reserva de memoria */
  A = (double*)calloc(N*n,sizeof(double));
  v = (double*)calloc(n+b,sizeof(double)); //tamaño solo aumenta por b
  w = (double*)calloc(n+b,sizeof(double)); //tamaño solo aumenta por b

  /* Inicializar datos */
  for(i=0; i<n; i++) A[i*N+(i+p)] = 2*b;
  for(i=0; i<n; i++) {
    for(j=0; j<N; j++) {
      if (i+p<j && abs(i+p-j)<=b) A[i*N+j] = -1.0;
    }
  }
  for(i=0; i<n; i++) v[i] = 1.0;

  /* Multiplicación de matrices */
  
  double t1,t2;
  
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  parmatvec(N,n,p,b,A,v,w);
  
  MPI_Barrier(MPI_COMM_WORLD);
  t2 = MPI_Wtime();  
  
  if (me==0) printf("Tiempo transcurrido: %f s.\n", t2-t1);
  
  /* Imprimir solución */
  for(i=0; i<n; i++) printf("w[%d] = %g\n", i+p, w[i]); //Solo imprimimos w sin la parte de arriba.

  free(A);
  free(v);
  free(w);

  MPI_Finalize();
  return 0;
}

