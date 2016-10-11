#include <malloc.h>
#include "array.h"

double*  Make1DArray(int m)
{
	double *a;
	a=(double*)malloc(sizeof(double)*m);
	return a;
}

double** Make2DArray(int m,int n)
{
	int i;
	double **a;
	a= (double**)malloc(sizeof(double*)*m);
	for(i=0;i<m;i++)
		a[i]= (double*)malloc(sizeof(double)*n);
	return a;
}

void free2DArray(double **a, int m)
{
	int i;
	for( i=0; i<m; i++ )
		free(a[i]);
	free(a);
}

