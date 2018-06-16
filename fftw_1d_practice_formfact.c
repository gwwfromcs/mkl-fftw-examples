// Compile: 
// On ICES machine, 
// 1. module load fftw3/3.3.4 intel/17.0; 
// 2. icc -I$FFTW_INC -L$FFTW_LIB -lfftw3 fftw_1d_practice_formfact.c -o fftw_1d_practice_formfact.exe
// 3. ./fftw_1d_practice.exe > data.txt
#include <complex.h>
#include <fftw3.h>
#include <math.h>

double gaussian(double mean, double stddev, double x) 
{ 
    double variance2 = stddev*stddev*2.0; 
    double term = x-mean; 
    double pi = 3.1415926535897932384626433832795;
    return exp(-(term*term)/variance2)/sqrt(pi*variance2); 
}  

int main(int argc, char **argv)
{
    const double pi = 3.1415926535897932384626433832795; 
    const double L = 1.0, mean = 0.5;
    int n2;
    int N = 50, S, M;
    float sigma = 0.01, R = 10.0, dx = 0.0005;
    fftw_plan p;
    int i, j;
    double or, oi, x, G;
    double complex sum;
    fftw_complex *in, *in2, *inff, \
      *out, *outff;
    double *box;
    
    if (argc == 5) {
       sscanf( *(argv+1), "%d", &N);
       sscanf( *(argv+2), "%d", &S);
       sscanf( *(argv+3), "%f", &R);
       sscanf( *(argv+4), "%f", &dx);
       if ( N < 5 || N > 1000 ) {
          printf (" N should between [5,1000]\n");
          return 0;
       }
       if ( S < 2 || S > 30 ) {
          printf (" S should between [2,300]\n");
          return 0;
       }
    }
    else {
       printf ("Usage: ./fftw_1d_practice_formfact N S R dx\n");
       printf (" N is the size of the grid. \n");
       return 0;
    }

    n2 = N/2;
    M = (int)(R/dx);
    in    = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
    in2   = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
    inff  = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
    out   = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
    outff = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex)); 
    box   = (double*) malloc(M * sizeof(double));

    printf ("# N  = %14d\n", N);
    printf ("# S  = %14d\n", S);
    printf ("# R  = %14.4f\n", R);
    printf ("# dx = %14.8f\n", dx);
    printf ("# M  = %14d\n", M);

    for(i = 0; i < N; i++) {
       x = (L*i)/(double)N;
       in[i] = 0.0;
       for (j = -S; j<=S; j++){
         in[i] = in[i] + gaussian(mean+j, sigma, x);
       }
       in2[i]   = 0.0;
       inff[i]  = 0.0;
       out[i]   = 0.0;
       outff[i] = 0.0;
    };

    p = fftw_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

    fftw_execute(p);
    
    for(i = 0; i < N; i++) {
       if ( i <= n2 ) {
           G = i*2.0*pi/L;
       }
       else {
           G = (i-N)*2.0*pi/L;
       }
       // printf ("%12.6f \n", G);
       sum = 0.0;
       for(j = 0; j < M; j++) {
           x = (j-M/2) * dx;
           box[j] = gaussian(0.0, sigma, x) * ( cos(G*x) + sin(G*x)*I );
       }
       for(j = 1; j < M/2; j++) {
           sum = sum + dx/3.0 * ( box[2*j-2] + box[2*j-1] * 4.0 + box[2*j] );
       }
       outff[i] = sum * ( cos(G*mean) + sin(G*mean)*I ) * (double)N;
    };
     
    p = fftw_plan_dft_1d( N, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE );
     
    fftw_execute(p);
    
    p = fftw_plan_dft_1d( N, outff, inff, FFTW_BACKWARD, FFTW_ESTIMATE );

    fftw_execute(p);

    /*
    printf("#    x   real(outff)   imag(outff)   real(out)   imag(out) \n");
    for(i = 0; i<N; i++) {
       or=creal(out[i]);
       oi=cimag(out[i]);
       if ( i <= n2 ) {
          printf("  %6d %19.12e %19.12e %19.12e %19.12e  \n", i, creal(outff[i]),cimag(outff[i]),or,oi);
       }
       else {
          printf("  %6d %19.12e %19.12e %19.12e %19.12e  \n", (i-N), creal(outff[i]),cimag(outff[i]),or,oi);
       }
    }
    */
    
    printf("#    x   real(in)  imag(in)   real(in2)   imag(in2)   real(inff)   imag(inff)\n");
    for(i = 0; i<N; i++) {
        printf("  %16.8f %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e  \n", (double)i/(double)N, \
          creal(in[i]), cimag(in[i]), creal(in2[i])/N, cimag(in2[i])/N, creal(inff[i])/N, cimag(inff[i])/N);
    }

    fftw_destroy_plan(p);
    return 0;
}
