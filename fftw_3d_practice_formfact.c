// Compilation on ICES:
// icc -I$FFTW_INC -L$FFTW_LIB -lfftw3 fftw_3d_practice_formfact.c -o fftw_3d_practice_formfact.exe
// This example shows how to use dynamically allocated 3d array for fftw
//
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

double gaussian(double center, double stddev, double x) 
{ 
    double variance2 = stddev*stddev*2.0; 
    double term = x-center; 
    // double pi = 3.1415926535897932384626433832795;
    const double pi = 3.1415926535897932384626433832795028841971693993751;
    return exp(-(term*term)/variance2)/sqrt(pi*variance2); 
}

double J0(double x)
{
    double small = 1.0e-8;
    if (x < small) {
        return 1.0;
    }
    else {
        return sin(x)/x;
    }
}

int main(int argc, char **argv)
{
    fftw_plan p;
    // const double pi = 3.1415926535897932384626433832795;
    const double pi = 3.1415926535897932384626433832795028841971693993751;
    const int L = 10;
    int i, j, k, N, M, Z, ngl, igl, ng, exist, mill[3];
    double gg, dx, Omega, theta, small = 1.0e-10, cutoff;
    double x,y,z,dist,r,sigma = 0.03, \
      sum, sum2, sumff, sumff_cut;
    double center[3], g[3];
    double complex *strf;
    double *gl, *formf, *box;
    int *gtogl;
    fftw_complex *in, *out, *outff, *outff_cut, *in2, \
      *inff, *inff_cut;
    
    // default value
    dx = 0.00001;
    Z  = 200000;
    center[0] = 0.5; center[1] = 0.5; center[2] = 0.5*L;

    // test bessel function J0(x)
    /*
    for (i = 0; i < Z; i++) {
       x = i * dx;
       printf (" %14.7e  %14.7e  %14.7e \n ",x, J0(x), gaussian(0.0, sigma, x)*x*x*J0(2*pi*5.0*x));
    }
    exit(0);
    */

    if (argc == 2) {
       sscanf(*(argv+1), " %d", &N);
       if ( N < 5 || N > 300 ) {
         printf (" N should between [5,300]\n");
         return 0;
       }
    }
    else {
       printf ("Usage: ./fftw_3d_practice_formfact N \n");
       printf (" N is the size of the grid in x and y direction. \n");
       return 0;
    }
    M = L * N;
    printf ("# N = %d M = %d Z = %d dx = %11.4f \n", N, M, Z, dx);
    in    = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    out   = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    outff = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    outff_cut \
          = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    inff  = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    inff_cut \
          = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    in2   = (fftw_complex *)fftw_malloc(N * N * M * sizeof(fftw_complex));
    for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
    for(k = 0; k < M; k++) {
       x = (double)i/N;    // x change from 0 to (N-1)/N
       y = (double)j/N;    // y change from 0 to (N-1)/N
       z = (double)k/M*L;  // z change from 0 to (M-1)/M*L
       dist = sqrt( (x-center[0]) * (x-center[0]) + \
                    (y-center[1]) * (y-center[1]) + \
                    (z-center[2]) * (z-center[2]) ); // distance from the grid point to center[3]
       in[k+M*(j+N*i)] = gaussian( 0.0, sigma, dist );
    }
     
    p = fftw_plan_dft_3d(N, N, M, \
            in, out, FFTW_FORWARD, FFTW_ESTIMATE);
     
    fftw_execute(p);
    // shells of g-vectors
    gl = (double *) malloc(N * N * M * sizeof(double));
    for (i = 0; i < N*N*M; i++) { gl[i] = 0.0; }
    gtogl = (int *) malloc(N * N * M * sizeof(int));
    for (i = 0; i < N*N*M; i++) { gtogl[i] = 0; }
    strf = (double complex *) malloc(N * N * M * sizeof(double complex));
    formf = (double *) malloc(N * N * M * sizeof(double));
    ngl = 0;
    for(i = 0; i < N; i++) {
      mill[0] = (i<=N/2)?(i):(i-N);
      g[0] = mill[0] * 2.0*pi/1.0;
      for(j = 0; j < N; j++) {
        mill[1] = (j<=N/2)?(j):(j-N);
        g[1] = mill[1] * 2.0*pi/1.0;
        for(k = 0; k < M; k++) {
           mill[2] = (k<=M/2)?(k):(k-M);
           g[2] = mill[2] * 2.0*pi/L;
           gg = sqrt ( g[0]*g[0] + g[1]*g[1] + g[2]*g[2] );
           ng = k+M*(j+N*i);
           exist = 0;
           for(igl = 0; igl < ngl; igl++) {
              if ( fabs(gg - gl[igl]) < small ) {
                 gtogl[ng] = igl;
                 exist = 1;
                 break;
              }
           }
           if ( !exist ) {
              ngl = ngl + 1;
              gtogl[ng] = ngl-1;
              gl[ngl-1] = gg;
           }
           // printf (" %5d %5d %5d %10.4e %10.4e %10.4e %10.4e %5d \n", mill[0], mill[1], mill[2], g[0], g[1], g[2], gg, gtogl[ng]);
           theta = g[0] * center[0] + \
             g[1] * center[1] + \
             g[2] * center[2];
           strf[ng] = cos(theta) + sin(theta)*I;

        } // k loop
      } // j loop
    } // i loop
    printf ( " # shells of gvecs: ngl = %9d \n", ngl );
    /*
    for (igl = 0;igl < ngl; igl++) {
        printf ( " %6d  %12.4e \n", igl, gl[igl] );
    }
    */
    Omega = 1.0 * 1.0 * L; // Volume of unit cell
    box = (double * ) malloc(Z * sizeof(double));
    for(igl = 0; igl < ngl; igl++) {
       if (gl[igl] <= small ) {
          for (j = 0; j < Z; j++) {
              x = j * dx;
              box[j] = gaussian(0.0, sigma, x) * x * x;
              // if ( j < 10000 && igl < 2) { printf ("%9.5e %9.5e \n", x, box[j]); }
          }
       }
       else {
          for (j = 0; j < Z; j++) {
              x = j * dx;
              box[j] = gaussian(0.0, sigma, x) * x * x * J0(gl[igl]*x);
              // if ( j < 10000 && igl < 2 ) { printf ("%9.5e %9.5e \n", x, box[j]); }
          }
       }
       sum = 0.0;
       // Integration with the Simpson's rule
       for (j = 1; j < Z/2; j++) {
           sum = sum + dx/3.0 * (box[2*j-2] + box[2*j-1]*4.0 + box[2*j]);
       }
       formf[igl] = sum * 4.0 * pi / Omega;
       // printf (" %9d %9.5f \n", igl, formf[igl]);
    }
    // printf ("    formf    strf \n");
    cutoff = (N/2 + 0.01)*2*pi/1.0;
    for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
    for (k = 0; k < M; k++) {
       mill[0] = (i<=N/2)?(i):(i-N);
       g[0] = mill[0] * 2.0*pi/1.0;
       mill[1] = (j<=N/2)?(j):(j-N);
       g[1] = mill[1] * 2.0*pi/1.0;
       mill[2] = (k<=M/2)?(k):(k-M);
       g[2] = mill[2] * 2.0*pi/L;
       gg = sqrt ( g[0]*g[0] + g[1]*g[1] + g[2]*g[2] );

       ng = k + M*(j+N*i);
       outff[ ng ] = formf[gtogl[ng]] * strf[ng] * N * N * M;
       if ( gg > cutoff ) { 
          outff_cut [ ng ] = 0.0 ;
       }
       else {
          outff_cut [ ng ] = outff [ ng ];
       }
       // printf( " %9.5e %9.5e \n", formf[gtogl[ng]], strf[ng]);
    }
    // compare out and outff
    /*
    printf ("#    out     outff \n");
    for(i = 0; i<N; i++)
    for(j = 0; j<N; j++)
    for(k = 0; k<M; k++) {
        mill[0] = (i<=N/2)?(i):(i-N);
        mill[1] = (j<=N/2)?(j):(j-N);
        mill[2] = (k<=M/2)?(k):(k-M);
        printf( " %9d %9d %9d %12.5e %12.5e %12.5e %12.5e \n", \
          mill[0], mill[1], mill[2], \
          creal( out[ k+M*(j+N*i) ] ), cimag( out[ k+M*(j+N*i) ] ), \
          creal( outff[ k+M*(j+N*i) ] ), cimag( outff[ k+M*(j+N*i) ] ) );
    }
    */
    p = fftw_plan_dft_3d(N, N, M, \
            out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    p = fftw_plan_dft_3d(N, N, M, \
            outff, inff, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    p = fftw_plan_dft_3d(N, N, M, \
        outff_cut, inff_cut, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    
    for(i = 0; i<N; i++)
    for(j = 0; j<N; j++)
    for(k = 0; k<M; k++) {
        x = (double)i/N; // change from 0 to 1
        y = (double)j/N; // change from 0 to 1
        z = (double)k/M*L; // change from 0 to L
        in2[ k+M*(j+N*i) ] = in2[ k+M*(j+N*i) ] / ((double)(N*N*M));
        inff[ k+M*(j+N*i) ] = inff[ k+M*(j+N*i) ] / ((double)(N*N*M));
        inff_cut[ k+M*(j+N*i) ] = inff_cut[ k+M*(j+N*i) ] / \
          ((double)(N*N*M));
        printf( "  %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e \n", \
          x, y, z, \
          creal( in[k+M*(j+N*i)] ), cimag( in[k+M*(j+N*i)] ),  \
          creal( in2[k+M*(j+N*i)] ), cimag( in2[k+M*(j+N*i)] ), \
          creal( inff[k+M*(j+N*i)] ), cimag( inff[k+M*(j+N*i)] ), \
          creal( inff_cut[k+M*(j+N*i)] ), cimag( inff_cut[k+M*(j+N*i)] ) );
    }
    printf ( " ==== " );
    
    
    for(k = 0; k<M; k++) {
       sum   = 0.0;
       sum2  = 0.0;
       sumff = 0.0;
       sumff_cut = 0.0;
       for(i = 0; i<N; i++)
       for(j = 0; j<N; j++) {
          sum = sum + creal(in[ k+M*(j+N*i) ]);
          sum2 = sum2 + creal(in2[ k+M*(j+N*i) ]); // / ((double)(N*N*M));
          sumff = sumff + creal(inff[ k+M*(j+N*i) ]); // / ((double)(N*N*M));
          sumff_cut = sumff_cut + \
             creal(inff_cut[ k+M*(j+N*i) ]); // / ((double)(N*N*M));
       }
       sum /= (double)(N*N);
       sum2 /= (double)(N*N);
       sumff /= (double)(N*N);
       sumff_cut /= (double)(N*N);
       z = (double)k/M*L;
       printf ( " %13.5e  %13.5e  %13.5e  %13.5e  %13.5e \n", \
          z, sum, sum2, sumff, sumff_cut );
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_free(outff);
    fftw_free(inff);
    fftw_free(in2);
    free(gl);
    free(gtogl);
    free(strf);
    free(formf);
    free(box);
    return 0;
}
