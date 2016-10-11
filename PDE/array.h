#ifndef ARRAY_H
#define ARRAY_H

#include <malloc.h>

double*  Make1DArray(int m);
double** Make2DArray(int m,int n);
void free2DArray(double **a, int m);

#endif

