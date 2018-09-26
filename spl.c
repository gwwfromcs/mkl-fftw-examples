#include <stdio.h>

// Numerical Recipe in C, Page. 115
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]) {
  /*
  Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
  x1 < x2 <...< xN , and given values yp1 and ypn for the first derivative of the interpolating
  function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
  the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
  ypn are equal to 10^30 or larger, the routine is signaled to set the corresponding boundary
  condition for a natural spline, with zero second derivative on that boundary
  */
  // ignore x[0] and y[0]
    int i,k;
    double p,qn,sig,un,*u;

    u = (double *)malloc(n*sizeof(double));
    // ignore u[0] start from u[1] to u[u-1]
    if (yp1 > 0.99e30)
        y2[1] = u[1] = 0.0;
    else {
        y2[1] = -0.5;
        u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) { 
        /*This is the decomposition loop of the tridiagonal algorithm.
        y2 and u are used for temporary
        storage of the decomposed
        factors. */
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)  // The upper boundary condition is set either to be
        qn=un=0.0;      // “natural”
    else {              // or else to have a specified first derivative.
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)         // This is the backsubstitution loop of the tridiagonal
        y2[k]=y2[k]*y2[k+1]+u[k];// algorithm.
    free(u);
}

// Numerical Recipe in C
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
/*
Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
and given the array y2a[1..n], which is the output from spline above, and given a value of
x, this routine returns a cubic-spline interpolated value y.
*/
{
    int klo,khi,k;
    double h,b,a;
    klo=1; /*We will find the right place in the table by means of
    bisection. This is optimal if sequential calls to this
    routine are at random values of x. If sequential calls
    are in order, and closely spaced, one would do better
    to store previous values of klo and khi and test if
    they remain appropriate on the next call. */
    khi=n;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    } // klo and khi now bracket the input value of x.
    h=xa[khi]-xa[klo];
    if (h == 0.0) {
      printf ( " Bad xa input to routine splint" );
      exit(0);
    }
    a=(xa[khi]-x)/h;  // The xa’s must be distinct.
    b=(x-xa[klo])/h; // Cubic spline polynomial is now evaluated.
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

// Numerical Recipe in C
void splint_array(double xa[], double ya[], double y2a[], int n, int m, double *x, double *y)
/*
Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
and given the array y2a[1..n], which is the output from spline above, and given an array
x[1..m], this routine returns a cubic-spline interpolated value y[1..m].
Assume x[1] < x[2] < x[3] < ... < x[m]
*/
{
    int klo,khi,k,i;
    double h,b,a;
    // check if xa[1] <= x[1] <= x[2] <= ... <= x[m] <= xa[n]
    if ( x[1] < xa[1] || x[1] > x[2] ) { 
      printf ("x[1] = %lf\n xa[1] = %lf\n x[2] = %lf \n", x[1], xa[1], x[2]);
      printf ("xa[1]<=x[1]<=x[2] doesnt holds" );
      exit(0);
    }
    for ( i = 2; i < m; i++ ) {
      if ( x[i] > x[i+1] ) { 
         printf ( "x[%d]=%lf\n x[%d]=%lf\n",i,x[i],i+1,x[i+1] );
         printf ( "x[%d] <= x[%d] doesnt holds ",i,i+1 );
         exit(0);
      }
    }
    if ( x[m] > xa[n] ) {
      printf ( "x[%d]=%lf\n xa[%d]=%lf\n",m,x[m],n,xa[n] );
      printf ( "x[%d] <= xa[%d] doesnt holds",m,n );
      exit(0);
    }
    klo=1; 
    khi=n;
    for(i = 1; i <= m; i++) {  
      if (x[i] > xa[khi]) {
        klo = khi;
        khi = n;
      }
      while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x[i]) { khi = k; }
        else { klo = k; }
      }
      h=xa[khi]-xa[klo];
      if (h == 0.0) {
        printf ( " Bad xa input to routine splint" );
        exit(0);
      }
      a=(xa[khi]-x[i])/h;  // The xa’s must be distinct.
      b=(x[i]-xa[klo])/h;  // Cubic spline polynomial is now evaluated.
      y[i]=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
    }
}

/*
 * input: func[1...n], we ignore func[0]
 *  rab[1...n]
 * output:
 *  asum
 */
double simpson(int n, double *func, double *rab)
{
    double f1, f2, f3, r12, asum;
    int i;
    asum = 0.0;
    r12 = 1.0/3.0;
    f3 = func[1] * rab[1] * r12;
    for ( i = 2; i <= n-1; i=i+2 ) {
        f1 = f3;
        f2 = func[i] * rab[i] * r12;
        f3 = func[i+1] * rab[i+1] * r12;
        asum = asum + f1 + 4.0 * f2 + f3;
    }
    return asum;
}

/*
 * input: func[1...n]
 *  dx
 * output: 
 *  asum
 */
double simpson2(int n, double *func, double dx)
{
    double f1, f2, f3, r12, asum;
    int i;
    asum = 0.0;
    r12 = 1.0/3.0;
    f3 = func[1] * dx * r12;
    for ( i = 2; i <= n-1; i=i+2 ) {
        f1 = f3;
        f2 = func[i] * dx * r12;
        f3 = func[i+1] * dx * r12;
        asum = asum + f1 + 4.0 * f2 + f3;
    }
    return asum;
}
