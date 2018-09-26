// declarations of auxilary functions 

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

void splint_array(double xa[], double ya[], double y2a[], int n, int m, double *x, double *y);

double simpson(int n, double *func, double *rab);

double simpson2(int n, double *func, double dx);
