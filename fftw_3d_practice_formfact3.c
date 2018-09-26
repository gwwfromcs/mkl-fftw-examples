// Compilation on ICES:
// icc -I$FFTW_INC -L$FFTW_LIB -lfftw3 spl.c fftw_3d_practice_formfact3.c -o fftw_3d_practice_formfact3.exe
// This example shows how to use dynamically allocated 3d array for fftw
//
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include "spl.h"

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

//     _____________ 
//    /|           /|
//   /            / |
//  /__|_________/  |
//  |           |   |
//  |  |        |   |
//  |           |   | 
//  |  |        |   |
//  |      @    |   | b
//  |  |        |   |
//  |           |   |
//  |  |        |   |
//  |           |   |
//  |  |_ _ _ _ |_ _|
//  |  /        |  /
//  |           | /  a
//  |/__________|/
//        a 
//

int main(int argc, char **argv)
{
    fftw_plan p;
    // const double pi = 3.1415926535897932384626433832795;
    const double pi = 3.1415926535897932384626433832795028841971693993751;
    const double twopi = 2.0*pi, fourpi = 4.0*pi;
    const double a = 4.5;   // in Bohr
    const int L = 10; 
    const double b = a * L;
    int i, j, k, N, M, Z, ngl, igl, ng, exist, mill[3], ncr;
    double gg, dx, Omega, theta, small = 1.0e-10, cutoff, tmpd1, tmpd2, tmpd3, tmpd4;
    double x,y,z,dist,r,sigma = 0.25, \
      sum, sum2, sumff, sumff_cut;
    double center[3], g[3];
    double complex *strf;
    double *gl, *formf, *integralbox, *box, *box2, *rhoc, *rr, *rhoc2, *rrab, *xx, \
      *rr_intp, *rhoc_intp, *rhoc2_intp;
    int *gtogl;
    fftw_complex *in, *out, *outff, *outff_cut, *in2, \
      *inff, *inff_cut;
    FILE *rhoinp;
    char filebuffer[100];
    
    // test bessel function J0(x)
    /*
    for (i = 0; i < Z; i++) {
       x = i * dx;
       printf (" %14.7e  %14.7e  %14.7e \n ",x, J0(x), gaussian(0.0, sigma, x)*x*x*J0(2*pi*5.0*x));
    }
    exit(0);
    */

    // default value
    dx = 0.00005; // in bohr, step size for integration
    Z  = 72000;  // number of points for integration
    center[0] = 0.5 * a; center[1] = 0.5 * a; center[2] = 0.5 * b; // in bohr
     
    // read command line argument and allocate space for FFT
    if (argc == 3) {
       sscanf(*(argv+1), " %d", &N);
       sscanf(*(argv+2), " %d", &M);
       if ( N < 5 || N > 300 ) {
         printf (" N should between [5,300]\n");
         return 0;
       }
    }
    else {
       printf ("Usage: ./fftw_3d_practice_formfact N M \n");
       printf (" N is the size of the grid in x and y direction. \n");
       printf (" M is the size of the grid in z direction. \n");
       return 0;
    }
     
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

    // read rhocore_r.dat
    ncr = 0;
    rhoinp = fopen("rhocore_r.dat", "r");
    // find out how many lines in rhocor_r.dat.
    while (( fgets(filebuffer, sizeof(filebuffer), rhoinp)) != NULL ) {
       ncr = ncr + 1;
    }
    // reposition rhoinp to the beginning of the file
    rewind(rhoinp);
    rhoc  = (double *) malloc ((ncr+1) * sizeof(double));
    rr    = (double *) malloc ((ncr+1) * sizeof(double));
    rrab  = (double *) malloc ((ncr+1) * sizeof(double));
    // I added rhocore(r=0) to rr_intp, rhoc_intp, rhoc2_intp
    rr_intp = (double *) malloc ((ncr+2) * sizeof(double));
    rhoc_intp = (double *) malloc ((ncr+2) * sizeof(double));
    rhoc2_intp = (double *) malloc ((ncr+2) * sizeof(double));
    for ( i = 1; i < ncr + 1; i++) {
      fscanf ( rhoinp, "%d %lf %lf %lf", &j, &tmpd1, &tmpd2, &tmpd3 );
      rr[i] = tmpd1; 
      rr_intp[i+1] = tmpd1; 
      rrab[i] = tmpd2; 
      rhoc[i] = tmpd3;
      rhoc_intp[i+1] = tmpd3;
      // printf ( " %6d  %24.12e %24.12e %24.12e \n ", j, rr[i], rrab[i], rhoc[i] );
    }
    // rhocore_r.dat doesn't have points at r=0, set rhoc[r=0] = rhoc[r(1)]
    rr_intp[1] = -small; rhoc_intp[1] = rhoc_intp[2]; 
    tmpd1 = 0.0; tmpd2 = 0.0;
    // setup the second derivative "rhoc2_intp" for splin interpolation
    spline(rr_intp, rhoc_intp, ncr+1, tmpd1, tmpd2, rhoc2_intp);
    // xx[1...Z], ignore xx[0]
    xx  = (double *) malloc ( (Z+1) * sizeof(double) );
    for ( i = 0; i < Z; i++) xx[i+1] = (double)i * dx;
    integralbox = (double *) malloc ( (Z+1) * sizeof(double) );
    // integralbox[1..Z] = rhoc[xx[1:Z]]
    splint_array(rr_intp, rhoc_intp, rhoc2_intp, ncr+1, Z, xx, integralbox);
    /* // test splin_array()
    for ( i = 0; i < Z; i++) {
       x = (double)i * dx;
       if ( x < rr_intp[1] ) {
          tmpd1 = rhoc[1];
       }
       else {
          splint(rr_intp, rhoc_intp, rhoc2_intp, ncr+1, x, &tmpd1);
       }
       printf ("%20.10e  %20.10e  %20.10e \n", x, tmpd1, integralbox[i+1]);
    }
    */
    // calculate rhoc (r) using realspace method 
    for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
    for(k = 0; k < M; k++) {
       x = (double)i/N * a;  // x change from 0 to (N-1)/N*a
       y = (double)j/N * a;  // y change from 0 to (N-1)/N*a
       z = (double)k/M * b;  // z change from 0 to (M-1)/M*b
       dist = sqrt( (x-center[0]) * (x-center[0]) + \
                    (y-center[1]) * (y-center[1]) + \
                    (z-center[2]) * (z-center[2]) ); // distance from the grid point to center[3]
       if ( dist < rr_intp[1] ) {
           tmpd1 = rhoc[1];
       } else {
           splint(rr_intp, rhoc_intp, rhoc2_intp, ncr+1, dist, &tmpd1);
       }
       in[k+M*(j+N*i)] = tmpd1;
    }
    p = fftw_plan_dft_3d(N, N, M, \
           in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Next, use fourier transform to calculate rhoc(G)
    // First find out how many shells of G-vectors 
    // gl stores the g-shells
    gl = (double *) malloc(N * N * M * sizeof(double));
    for (i = 0; i < N*N*M; i++) { gl[i] = 0.0; }
    // gtogl maps a g-vector to its index of gl
    gtogl = (int *) malloc(N * N * M * sizeof(int));
    for (i = 0; i < N*N*M; i++) { gtogl[i] = 0; }
    // structure factor
    strf = (double complex *) malloc(N * N * M * sizeof(double complex));
    // form factor
    formf = (double *) malloc(N * N * M * sizeof(double));
    ngl = 0;
    for (i = 0; i < N; i++) {
      mill[0] = (i<=N/2)?(i):(i-N);
      g[0] = mill[0] * twopi/a;
      for(j = 0; j < N; j++) {
        mill[1] = (j<=N/2)?(j):(j-N);
        g[1] = mill[1] * twopi/a;
        for(k = 0; k < M; k++) {
           mill[2] = (k<=M/2)?(k):(k-M);
           g[2] = mill[2] * twopi/b;
           gg = sqrt ( g[0]*g[0] + g[1]*g[1] + g[2]*g[2] ); // In Rydberg
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
    /* // print out hte shells of g-vectors
    for (igl = 0;igl < ngl; igl++) {
        printf ( " %6d  %12.4e \n", igl, gl[igl] );
    }
    */
    Omega = a*a*b; // Volume of unit cell
    box = (double *) malloc ( (ncr+1)*sizeof(double));
    box2 = (double *) malloc ( (Z+1)*sizeof(double));
    for(igl = 0; igl < ngl; igl++) {
        // integrate with finer grid
       if (gl[igl] <= small ) {
          for (j = 1; j <= Z; j++) {
              x = xx[j];
              box2[j] = integralbox[j] * x * x;
          }
       }
       else {
          for (j = 1; j <= Z; j++) {
              x = xx[j];
              box2[j] = integralbox[j] * x * x * J0(gl[igl]*x);
          }
       }
       //sum = simpson2(Z, box2, dx);
       sum = 0.0;
       for ( j = 1; j < Z/2; j++) {
           sum = sum + dx/3.0 * (box2[2*j-2+1] + box2[2*j-1+1]*4.0 + box2[2*j+1]);
       }

       // printf ("%20.12e  %20.12e \n", gl[igl], sum);
       /*
       // integrate as QE
       if (gl[igl] <= small ) {
          for (j = 1; j <= ncr; j++) {
              x = rr[j];
              box[j] = rhoc[j] * x * x;
              // printf ("%7d %15.8f %20.12e %20.12e\n",j, x, box[j], rhoc[j]);
          }
       }
       else {
          for (j = 1; j <= ncr; j++) {
              x = rr[j];
              box[j] = rhoc[j] * x * x * J0(gl[igl]*x);
              // printf ("%7d %15.8f %20.12e %20.12e\n",j, x, box[j], rhoc[j]);
          }
       }
       sum = simpson(ncr, box, rrab);
       printf (" %20.12e\n" sum);
       */
       formf[igl] = sum * fourpi / Omega; // divided by Omega
       // print out four factor
       // If there is only one atom in the unit cell, then |formf(|G|)| = V(G)
       printf (" %8d %15.8e %15.8e \n", igl, gl[igl]*gl[igl], formf[igl] );
    }
    // printf ("    formf    strf \n");
    exit(0);
    cutoff = (N/2 + 0.01)*twopi/a;
    for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
    for (k = 0; k < M; k++) {
       mill[0] = (i<=N/2)?(i):(i-N);
       g[0] = mill[0] * twopi/a;
       mill[1] = (j<=N/2)?(j):(j-N);
       g[1] = mill[1] * twopi/a;
       mill[2] = (k<=M/2)?(k):(k-M);
       g[2] = mill[2] * twopi/b;
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
       // printf( " %10d %10d %10d %17.8e %17.8e \n", mill[0], mill[1], mill[2], outff[ng]/N/N/M, formf[gtogl[ng]] );
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
        x = (double)i/N*a; // change from 0 to a
        y = (double)j/N*a; // change from 0 to a
        z = (double)k/M*b; // change from 0 to b
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
    printf ( " ==== \n" );
    
    for(k = 0; k < M; k++) {
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
       z = (double)k/M*b;
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
    free(integralbox);
    free(rhoc);
    free(rr);
    free(rrab);
    free(rr_intp);
    free(rhoc_intp);
    free(rhoc2_intp);
    free(xx);
    free(box);
    free(box2);
        
    return 0;
}
