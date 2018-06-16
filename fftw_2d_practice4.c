// Compilation on ICES:
// icc -I$FFTW_INC -L$FFTW_LIB -lfftw3 fftw_2d_practice4.c -o fftw_2d_practice4.exe
// This example shows how to use dynamically allocated 2d array for fftw
//
#include <complex.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

int main(int argc, char **argv)
{
    fftw_plan p;
    const double pi = 3.141592653589793238462643;
    const int L = 10;
    int i, j, N, M;
    char zz;
    double radius;
    double x,y,r,sigma = 0.08;
    fftw_complex *in;
    fftw_complex *out;
    fftw_complex *in2;
    
    if (argc == 3) {
       sscanf(*(argv+1), " %d", &N);
       sscanf(*(argv+2), " %c", &zz);
       if ( N < 5 || N > 300 ) {
         printf (" N should between [5,300]\n");
         return 0;
       }
       if ( zz != 'y' && zz != 'n' ) {
         printf (" zz should be either y or n\n");
         return 0;
       }
    }
    else {
       printf ("Usage: ./fftw_2d_practice4 N z \n");
       printf (" where z = 'y' or 'n' \n");
       printf (" 'y' is for set some of FFT[in] to zero. \n");
       printf (" N is the size of the grid. \n");
       return 0;
    }
    M = L * N;
    printf ("# N = %d, M = %d \n", N, M);
    printf ("# zz = %c\n", zz);
    radius = (double)N/2+0.01;
    in  = (fftw_complex *)fftw_malloc(N*M* sizeof(fftw_complex));
    out = (fftw_complex *)fftw_malloc(N*M* sizeof(fftw_complex));
    in2 = (fftw_complex *)fftw_malloc(N*M* sizeof(fftw_complex));

    for(i = 0; i<N; i++)
    for(j = 0; j<M; j++){
       x = (double)i/N; // change from 0 to 1.0
       y = (double)j/N; // change from 0 to L
       // Accessing in[i][j]
       /*
       in[j+M*i] = sqrt(1./(2.*pi*sigma)) * exp((-pow(x,2)-pow(y,2))/2.0/pow(sigma,2)) +  \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x,2)-pow(y-1.0,2))/2.0/pow(sigma,2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x-1.0,2)-pow(y,2))/2.0/pow(sigma,2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x-1.0,2)-pow(y-1.0,2))/2.0/pow(sigma,2)) ;
       */
       /*
       in[j+M*i] = sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.50 * L, 2))/2.0/pow(sigma, 2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.55 * L, 2))/2.0/pow(sigma, 2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.60 * L, 2))/2.0/pow(sigma, 2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.65 * L, 2))/2.0/pow(sigma, 2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.70 * L, 2))/2.0/pow(sigma, 2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.75 * L, 2))/2.0/pow(sigma, 2)) ;
       */
       //
       in[j+M*i] = sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.50 * L, 2))/2.0/pow(sigma, 2)) ;
    }
     
    p = fftw_plan_dft_2d(N, M, \
            in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    // zero out some elements of out[i][j] 
    if ( zz == 'y' ) {
      /* cutoff radius
      for(i = 0; i<N; i++)
      for(j = 0; j<M; j++){
        if ( sqrt((double)(i-N)*(i-N)+(double)j*j/(L*L)) > radius &&
             sqrt((double)i*i+(double)(j-M)*(j-M)/(L*L)) > radius &&
             sqrt((double)(i-N)*(i-N)+(double)(j-M)*(j-M)/(L*L)) > radius &&
             sqrt((double)i*i+(double)j*j/(L*L)) > radius ) {
          out[j+M*i] = 0.0;
        }
      }
         cutoff radius */
      /* middle planes */
      for(i = 0; i<N; i++)
      for(j = 0; j<M; j++){
        if ( i == N/2+1 ) {out[j+M*i] = 0.0;}
        if ( i%2 == 0 && i == N/2 ) {out[j+M*i] = 0.0;}
        if ( j == M/2+1 ) {out[j+M*i] = 0.0;}
        if ( j%2 == 0 && j == M/2 ) {out[j+M*i] = 0.0;}
      }
      /* middle planes */
    }

    p = fftw_plan_dft_2d(N, M, \
            out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    for(i = 0; i<N; i++)
    for(j = 0; j<M; j++){
        x = (double)i/N; // change from 0 to 1
        y = (double)j/N; // change from 0 to L
        in2[j+M*i] = in2[j+M*i]/((double)(N*M));
        printf("  %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e\n", \
          x, y, \
          creal(in[j+M*i]), cimag(in[j+M*i]), \
          creal(out[j+M*i]),cimag(out[j+M*i]), \
          creal(in2[j+M*i]),cimag(in2[j+M*i]));
    }


    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_free(in2);

    return 0;
}
