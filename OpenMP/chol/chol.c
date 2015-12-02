#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./quark-0.9.0/quark.h"

/* 
 * Convertir matriz cuadrada entre formato fortran estandar (por columnas,
 * lda=n) y formato de matriz tile con tamaño de bloque bs.
 * A es entrada, B es salida
 * Si tofrom=0, A es la de formato fortran, B es de tipo tile
 * Si tofrom=1, A es la de formato tile, B es de tipo fortran
 */
void convtile(double *A, double *B, int n, int bs, int tofrom)
{
  int i, j, ii, jj, i2, j2, nb=n/bs;

  for (i=0; i<nb; i++) {
    for (j=0; j<nb; j++) {
      for (jj=0; jj<bs; jj++) {
        if (tofrom) memcpy(B+i*bs+(j*bs+jj)*n, A+(i+j*nb)*bs*bs+jj*bs, bs*sizeof(double));
        else memcpy(B+(i+j*nb)*bs*bs+jj*bs, A+i*bs+(j*bs+jj)*n, bs*sizeof(double));
      }
    }
  }
}

/* 
 * Cholesky por bloques, almacenamiento fortran
 */
void choleskyblk(double *A, int n, int bs, int *info)
{
  int i, j, k, nb=n/bs, lda=n;
  const double done = 1.0, dmone = -1.0, dzero = 0.0;

  for (k=0; k<nb; k++) {
    dpotrf_("L", &bs, A+k*bs+k*bs*n, &lda, info, 1);
    if (*info != 0) return;
    for (j=k+1; j<nb; j++)
      dtrsm_("R", "L", "T", "N", &bs, &bs, &done, A+k*bs+k*bs*n, &lda, A+j*bs+k*bs*n, &lda, 1, 1, 1, 1);
    for (j=k+1; j<nb; j++) {
      for (i=j+1; i<nb; i++)
        dgemm_("N", "T", &bs, &bs, &bs, &dmone, A+i*bs+k*bs*n, &lda, A+j*bs+k*bs*n, &lda, &done, A+i*bs+j*bs*n, &lda, 1, 1);
      dsyrk_("L", "N", &bs, &bs, &dmone, A+j*bs+k*bs*n, &lda, &done, A+j*bs+j*bs*n, &lda);
    }
  }
}

/* ===================================================================================================================== */

/* Create a task wrapper to dpotrf, usable by the QUARK runtime.
 * Basically, you need to unpack the arguments from QUARK and call the
 * routine */
void dpotrf_quark_task( Quark *quark )
{
	printf("dpotrf_quark_task\n");
	
    double *A;
    int bs, lda;
    int * info;
    
	//printf("Unpack\n");
    quark_unpack_args_4( quark, bs, A, lda, info);
    
	//printf("dpotrf\n");
    dpotrf_("L", &bs, A, &lda, info, 1);
	//printf("POST dpotrf\n");
}


/* Create a call that will insert a dpotrf_ task into the QUARK
 * runtime.  Later, when dependencies are statisfied, the runtime will
 * execute this task.  The arguments to dpotrf_ are specified and
 * passed to QUARK here. */
void dpotrf_quark_call( Quark *quark, int bs, double *A, int lda, int * info)
{
	//printf("=> dpotrf_quark_call dpotrf\n");
    QUARK_Insert_Task( quark, dpotrf_quark_task, NULL,
	    				sizeof(int), &bs, VALUE,
						sizeof(double)*bs*bs, A, INOUT,
						sizeof(int), &lda, VALUE,
						sizeof(int), info, INOUT,						
						0 );
}

/* ===================================================================================================================== */

void dtrsm_quark_task( Quark *quark )
{
	printf("dtrsm_quark_task\n");
	
    double *A1;
	double *A2;
    int bs, lda;
    quark_unpack_args_4( quark, bs, A1, lda, A2);    
	const double done = 1.0;
	
	dtrsm_("R", "L", "T", "N", &bs, &bs, &done, A1, &lda, A2, &lda, 1, 1, 1, 1);    

}

void dtrsm_quark_call( Quark *quark, int bs, double *A1, int lda, double *A2)
{
    QUARK_Insert_Task( quark, dtrsm_quark_task, NULL,
	    				sizeof(int), &bs, VALUE,
						sizeof(double)*bs*bs, A1, INOUT,
						sizeof(int), &lda, VALUE,
						sizeof(double)*bs*bs, A2, INOUT,						
						0 );
}

/* ===================================================================================================================== */

void dgemm_quark_task( Quark *quark )
{
	printf("dgemm_quark_task\n");
	
    double *A1;
	double *A2;
	double *A3;
    int bs, lda;
    quark_unpack_args_5( quark, bs, A1, lda, A2, A3);    
	const double done = 1.0, dmone = -1.0, dzero = 0.0;
	
	dgemm_("N", "T", &bs, &bs, &bs, &dmone, A1, &lda, A2, &lda, &done, A3, &lda, 1, 1);
	
}

void dgemm_quark_call( Quark *quark, int bs, double *A1, int lda, double *A2, double *A3)
{
    QUARK_Insert_Task( quark, dgemm_quark_task, NULL,
	    				sizeof(int), &bs, VALUE,
						sizeof(double)*bs*bs, A1, INOUT,
						sizeof(int), &lda, VALUE,
						sizeof(double)*bs*bs, A2, INOUT,
						sizeof(double)*bs*bs, A3, INOUT,
						0 );
}

/* ===================================================================================================================== */

void dsyrk_quark_task( Quark *quark )
{
	printf("dgemm_quark_task\n");
	
    double *A1;
	double *A2;
    int bs, lda;
    quark_unpack_args_4( quark, bs, A1, lda, A2);    
	const double done = 1.0, dmone = -1.0, dzero = 0.0;
	
	dsyrk_("L", "N", &bs, &bs, &dmone, A1, &lda, &done, A2, &lda);
	
}

void dsyrk_quark_call( Quark *quark, int bs, double *A1, int lda, double *A2)
{
    QUARK_Insert_Task( quark, dsyrk_quark_task, NULL,
	    				sizeof(int), &bs, VALUE,
						sizeof(double)*bs*bs, A1, INOUT,
						sizeof(int), &lda, VALUE,
						sizeof(double)*bs*bs, A2, INOUT,
						0 );
}


/* ===================================================================================================================== */

/* 
 * Cholesky por bloques, almacenamiento tile
 */
void choleskytile(double *A, int n, int bs, int *info, Quark *quark)
{
  int i, j, k, nb=n/bs, lda=bs, bs2=bs*bs;
  const double done = 1.0, dmone = -1.0, dzero = 0.0;

  for (k=0; k<nb; k++) {
	
	dpotrf_quark_call(quark, bs,A+(k+k*nb)*bs2,lda,info);
    //dpotrf_("L", &bs, A+(k+k*nb)*bs2, &lda, info, 1);
	
	//QUARK_Barrier( quark );
    
    if (*info != 0) return;
    
    for (j=k+1; j<nb; j++){
		dtrsm_quark_call(quark, bs, A+(k+k*nb)*bs2, lda, A+(j+k*nb)*bs2);

		//dtrsm_("R", "L", "T", "N", &bs, &bs, &done, A+(k+k*nb)*bs2, &lda, A+(j+k*nb)*bs2, &lda, 1, 1, 1, 1);    
    }
    
	//QUARK_Barrier( quark );

    
    for (j=k+1; j<nb; j++) {
		for (i=j+1; i<nb; i++){  
		  dgemm_quark_call( quark, bs, A+(i+k*nb)*bs2, lda, A+(j+k*nb)*bs2, A+(i+j*nb)*bs2);
		
		  //dgemm_("N", "T", &bs, &bs, &bs, &dmone, A+(i+k*nb)*bs2, &lda, A+(j+k*nb)*bs2, &lda, &done, A+(i+j*nb)*bs2, &lda, 1, 1);    
		
		  //QUARK_Barrier( quark );
		}
		
		dsyrk_quark_call( quark, bs, A+(j+k*nb)*bs2, lda, A+(j+j*nb)*bs2);
		
		//dsyrk_("L", "N", &bs, &bs, &dmone, A+(j+k*nb)*bs2, &lda, &done, A+(j+j*nb)*bs2, &lda);
		//QUARK_Barrier( quark );
    }
  }
}

int main(int argc, char **argv)
{
  int i, j, n=640, bs=64, ft=0, info, nthreads=2;
  double *A, *B, *x, *b, err;
  const double done = 1.0, dzero = 0.0;
  const int one = 1;
  
  struct timeval t0, t1;
  double time;

  /* Extracción de argumentos */
  if (argc >= 2) {
    if ((n = atoi(argv[1])) < 0) n = 640;
  }
  if (argc >= 3) {
    if ((bs = atoi(argv[2])) < 0) bs = 64;
  }
  if (n<bs || n%bs) {
    printf("ERROR: el tamaño de la matriz ha de ser divisible por bs\n");
    exit(1);
  }
  if (argc >= 4) {
    if ((nthreads = atoi(argv[3])) < 0) nthreads = 640;
  }
  if (argc >= 5) ft = atoi(argv[4]);  /* ft=0 tile, ft=1 fortran */ 
  printf("Cholesky de matriz de orden %d, tamaño de bloque=%d %s\n",n,bs,ft?"FORTRAN":"TILE");
  

  /* Creación de las matrices */
  A = (double*)malloc(n*n*sizeof(double));
  B = (double*)malloc(n*n*sizeof(double));
  x = (double*)malloc(n*sizeof(double));
  b = (double*)malloc(n*sizeof(double));

  //Quark start threads
  Quark *quark = QUARK_New(nthreads);


  /* Inicializar matriz aleatoria */
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      A[i+n*j] = drand48();
      A[j+n*i] = A[i+n*j];
    }
  }
  for(i=0; i<n; i++) A[i+i*n] += (double)n;  /* necesario para que sea definida */

  /* Calcular b=A*ones(n,1), para comprobar posteriormente */
  for(i=0; i<n; i++) x[i] = 1.0;
  dgemv_("N", &n, &n, &done, A, &n, x, &one, &dzero, b, &one, 1);


  /* Cholesky */
  if (ft) choleskyblk(A, n, bs, &info);
  else {
    convtile(A, B, n, bs, 0);

	//get time
	gettimeofday (&t0, NULL);

    choleskytile(B, n, bs, &info, quark);
	
	//get time
	gettimeofday (&t1, NULL);
    
    convtile(B, A, n, bs, 1);
  }
  
  //Quark eliminate threads
  QUARK_Delete(quark);

  if (info != 0) {
    printf("ERROR: falló la factorización de Cholesky\n");
    exit(1);
  }

  /* Comprobar solución */
  dpotrs_("L", &n, &one, A, &n, b, &n, &info, 1);
  if (info != 0) {
    printf("ERROR: falló la resolución triangular\n");
    exit(1);
  }
  err = 0.0;
  for (i=0; i<n; i++) err += fabs(1.0-b[i]);
  printf("Error=%g\n",err);

  time = (t1.tv_sec-t0.tv_sec)+(t1.tv_usec-t0.tv_usec)/1000000.0;

  printf("%d;%d;%d;%f\n",n,bs,nthreads,time);


  free(A); free(B); free(x); free(b);
  return 0;
}

