#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
void jacobi_poisson(int N,int M,double *x, double *b)
{
  int i, j, k, conv, maxit=10000;
  double *t, s, tol=1e-6;

  t = (double*)calloc(N*M,sizeof(double));

  k = 0;
  conv = 0;

  while(!conv && k<maxit) {

    /* línea superior de la malla */
    i=0;
    t[i*M+0] = (b[i*M+0] + x[(i+1)*M+0] + x[i*M+(0+1)])/4.0;
    for(j=1; j<M-1; j++) {
      t[i*M+j] = (b[i*M+j] + x[(i+1)*M+j] + x[i*M+(j+1)] + x[i*M+(j-1)])/4.0;
    }
    t[i*M+(M-1)] = (b[i*M+(M-1)] + x[(i+1)*M+(M-1)] + x[i*M+((M-1)-1)])/4.0;

    /* interior de la malla */
    for(i=1; i<N-1; i++) {
      t[i*M+0] = (b[i*M+0] + x[(i+1)*M+0] + x[(i-1)*M+0] + x[i*M+(0+1)])/4.0;
      for(j=1; j<M-1; j++) {
        t[i*M+j] = (b[i*M+j] + x[(i+1)*M+j] + x[(i-1)*M+j] + x[i*M+(j+1)] + x[i*M+(j-1)])/4.0;
      }
      t[i*M+(M-1)] = (b[i*M+(M-1)] + x[(i+1)*M+(M-1)] + x[(i-1)*M+(M-1)] + x[i*M+((M-1)-1)])/4.0;
    }

    /* línea inferior de la malla */
    i=N-1;
    t[i*M+0] = (b[i*M+0] + x[(i-1)*M+0] + x[i*M+(0+1)])/4.0;
    for(j=1; j<M-1; j++) {
      t[i*M+j] = (b[i*M+j] + x[(i-1)*M+j] + x[i*M+(j+1)] + x[i*M+(j-1)])/4.0;
    }
    t[i*M+(M-1)] = (b[i*M+(M-1)] + x[(i-1)*M+(M-1)] + x[i*M+((M-1)-1)])/4.0;

    /* criterio de parada: ||x_{k}-x_{k+1}||<tol */
    s = 0.0;
    for(i=0; i<N; i++) {
      for(j=0; j<M; j++) {
        s += (x[i*M+j]-t[i*M+j])*(x[i*M+j]-t[i*M+j]);
      }
    }
    conv = (sqrt(s)<tol);
    printf("Error en iteración %d: %g\n", k, sqrt(s));

    /* siguiente iteración */
    k = k+1;
    for(i=0; i<N; i++) {
      for(j=0; j<M; j++) {
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

  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de N */
     if ((N = atoi(argv[1])) < 0) N = 50;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de M */
     if ((M = atoi(argv[2])) < 0) M = 1;
  }

  /* Reserva de memoria */
  x = (double*)calloc(N*M,sizeof(double));
  b = (double*)calloc(N*M,sizeof(double));

  /* Inicializar datos */
  for(i=0; i<N; i++) {
    for(j=0; j<M; j++) {
      b[i*M+j] = h*h*f;  /* suponemos que la función f es constante en todo el dominio */
    }
  }

  /* Resolución del sistema por el método de Jacobi */
  jacobi_poisson(N,M,x,b);
  
  /* Imprimir solución */
  for(i=0; i<N; i++) {
    for(j=0; j<M; j++) {
      printf("%g ", x[i*M+j]);
    }
    printf("\n");
  }

  free(x);
  free(b);

  return 0;
}

