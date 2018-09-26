
// Compile: 
// On ICES machine, 
// 1. module load fftw3/3.3.4 intel/17.0; 
// 2. icc -I$FFTW_INC -L$FFTW_LIB -lfftw3 fftw_1d_practice2.c -o fftw_1d_practice2.exe
// 3. ./fftw_1d_practice2.exe > data.txt
//
#include <complex.h>
#include <fftw3.h>
#include <math.h>

int main()
{
    fftw_plan p;
    const double pi = 3.141592653589793238462643;
    int i, N = 2000, nsquare = 10, offset=1;
    double sigma = 0.01, or, oi;
    fftw_complex in[N];
    fftw_complex out[N];
    fftw_complex in2[N];

    for(i = 0; i<N; i++) {
      // Gaussian function
      in[i] = sqrt(1./(2.*pi*sigma)) * exp(-((double)i/N-0.0)*((double)i/N-0.0)/(2.0*sigma*sigma)) +  \
        sqrt(1./(2.*pi*sigma)) * exp(-((double)i/N-1.0)*((double)i/N-1.0)/(2.0*sigma*sigma));
    };

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);
    
    p = fftw_plan_dft_1d(N, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    printf("#  out = fft(in) \n");
    printf("#    x   real(in(x))   imag(in(x))   real(out)   imag(out)  real(in2(x))  imag(in2(x)) \n");
    for(i = 0; i<N; i++) {
        in2[i] = in2[i]/((double)N);
        or=creal(out[i]);
        oi=cimag(out[i]);
        //printf( "  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n", \
        //  (double)i/N, creal(in[i]), cimag(in[i]), or, oi, creal(in2[i]), cimag(in2[i]) );
        printf( "  %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e \n", \
          (double)i/N, creal(in[i]), cimag(in[i]), or, oi, creal(in2[i]), cimag(in2[i]) );
    }

    fftw_destroy_plan(p);
    return 0;

}
