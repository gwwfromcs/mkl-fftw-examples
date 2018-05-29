#include <complex.h>
#include <fftw3.h>
#include <math.h>

int main()
{
    fftw_plan p;
    int i, j, N = 40;
    float x,y,r,sigma=0.05;
    fftw_complex in[N][N];
    fftw_complex out[N][N];

    for(i = 0; i<N; i++)
    for(j = 0; j<N; j++){
        x = (float)i/N;
        y = (float)j/N;
        in[i][j] = sqrt(1./(2.*3.1415926*sigma)) * exp((-pow(x,2)-pow(y,2))/2.0/pow(sigma,2)) +  \
                   sqrt(1./(2.*3.1415926*sigma)) * exp((-pow(x,2)-pow(y-1.0,2))/2.0/pow(sigma,2)) + \
                   sqrt(1./(2.*3.1415926*sigma)) * exp((-pow(x-1.0,2)-pow(y,2))/2.0/pow(sigma,2)) + \
                   sqrt(1./(2.*3.1415926*sigma)) * exp((-pow(x-1.0,2)-pow(y-1.0,2))/2.0/pow(sigma,2)) ;
    }

    p = fftw_plan_dft_2d(N, N, \
            &in[0][0], &out[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    for(i = 0; i<N; i++)
    for(j = 0; j<N; j++){
        x = (float)i/N;
        y = (float)j/N;
        printf("  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", x,y,creal(in[i][j]),cimag(in[i][j]),creal(out[i][j]),cimag(out[i][j]));
    }
    
    fftw_destroy_plan(p);
    return 0;
}
