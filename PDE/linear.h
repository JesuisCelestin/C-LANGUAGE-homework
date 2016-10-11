#ifndef LINEAR_H
#define LINEAR_H

// solve equations Ax=b, A is n*n matrix, note: A is a dynamically allocated 2d arrays
int GaussPivoting( double **A, double b[], int n, double x[] );

// Jacobi iteration
// int JacobiIter( double **A, double b[], double x[], int n );
// Gauss-Seidel iteration
int GaussSeidelIter( double **A, double b[], double x0[], int n, double x[] );

int chase(double a[],double b[],double c[], double f[], double yy[],int n);

#endif
