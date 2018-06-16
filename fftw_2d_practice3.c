// Compilation on ICES:
// icc -I$FFTW_INC -L$FFTW_LIB -lfftw3 fftw_2d_practice3.c -o fftw_2d_practice3.exe
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
    int i, j, N;
    char zz;
    double radius;
    double x,y,r,sigma = 0.05;
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
       printf ("Usage: ./fftw_2d_practice3 N z \n");
       printf (" where z = 'y' or 'n' \n");
       printf (" 'y' is for set some of FFT[in] to zero. \n");
       printf (" N is the size of the grid. \n");
       return 0;
    }

    printf ("# N = %d\n", N);
    printf ("# zz = %c\n", zz);
    radius = (double)N/2+0.01;
    in  = (fftw_complex *)fftw_malloc(N*N* sizeof(fftw_complex));
    out = (fftw_complex *)fftw_malloc(N*N* sizeof(fftw_complex));
    in2 = (fftw_complex *)fftw_malloc(N*N* sizeof(fftw_complex));

    for(i = 0; i<N; i++)
    for(j = 0; j<N; j++){
       x = (double)i/N;
       y = (double)j/N;
       // Accessing in[i][j]
       /*
       in[j+N*i] = sqrt(1./(2.*pi*sigma)) * exp((-pow(x,2)-pow(y,2))/2.0/pow(sigma,2)) +  \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x,2)-pow(y-1.0,2))/2.0/pow(sigma,2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x-1.0,2)-pow(y,2))/2.0/pow(sigma,2)) + \
                   sqrt(1./(2.*pi*sigma)) * exp((-pow(x-1.0,2)-pow(y-1.0,2))/2.0/pow(sigma,2)) ;
       */
       in[j+N*i] = sqrt(1./(2.*pi*sigma)) * exp((-pow(x - 0.5, 2)-pow(y - 0.5, 2))/2.0/pow(sigma, 2));
    }

    p = fftw_plan_dft_2d(N, N, \
            in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    // zero out some elements of out[i][j] 
    /*
    if ( zz == 'y' ) {
      for(i = 0; i<N; i++)
      for(j = 0; j<N; j++){
        if ( sqrt((double)(i-N)*(i-N)+(double)j*j) > radius &&
             sqrt((double)i*i+(double)(j-N)*(j-N)) > radius &&
             sqrt((double)(i-N)*(i-N)+(double)(j-N)*(j-N)) > radius &&
             sqrt((double)i*i+(double)j*j) > radius ) {
          out[j+N*i] = 0.0;
        }
      }
    }
    */

    p = fftw_plan_dft_2d(N, N, \
            out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    for(i = 0; i<N; i++)
    for(j = 0; j<N; j++){
        x = (double)i/N;
        y = (double)j/N;
        in2[j+N*i] = in2[j+N*i]/((double)(N*N));
        printf("  %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e\n", \
          x, y, \
          creal(in[j+N*i]), cimag(in[j+N*i]), \
          creal(out[j+N*i]),cimag(out[j+N*i]), \
          creal(in2[j+N*i]),cimag(in2[j+N*i]));
    }


    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_free(in2);

    return 0;
}
