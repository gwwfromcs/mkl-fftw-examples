/* C source code is found in dgemm_example.c */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include <time.h>

int main(int argc, char *argv[])
{
    double *A, *B, *C;
    int m, n, k, i, j;
    double alpha, beta;
	clock_t alloc_timer, cal_timer;
	int mtxsize;

	mtxsize = atoi(argv[1]);
	if (mtxsize <= 0)
	{
		printf("input argument is wrong");
		return 1;
	}

    m = mtxsize, k = mtxsize, n = mtxsize;
    alpha = 1.0; beta = 0.0;
    
	alloc_timer = clock();
    A = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
    B = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
    C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
    if (A == NULL || B == NULL || C == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(A);
      mkl_free(B);
      mkl_free(C);
      return 1;
    }

	srand(time(NULL));
    for (i = 0; i < (m*k); i++) {
        A[i] = ((double)(rand()-RAND_MAX/2.0))/((double)RAND_MAX);
    }

    for (i = 0; i < (k*n); i++) {
        B[i] = ((double)(rand()-RAND_MAX/2.0))/((double)RAND_MAX);
    }

    for (i = 0; i < (m*n); i++) {
        C[i] = 0.0;
    }
	alloc_timer = clock()-alloc_timer;

	cal_timer = clock();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A, k, B, n, beta, C, n);
	cal_timer = clock()-cal_timer;

	printf("Allocating time: %12.6f \n",(double)alloc_timer/CLOCKS_PER_SEC);  
	printf("Calculating time: %12.6f \n",(double)cal_timer/CLOCKS_PER_SEC);

    mkl_free(A);
    mkl_free(B);
    mkl_free(C);

    return 0;
}
