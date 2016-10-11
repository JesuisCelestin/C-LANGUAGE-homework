#include <math.h>
#include <stdio.h>
#include "linear.h"

// trace method in solving the tridiagonal equations
// [ b c ]
// [ a b ] yy = f
int chase(double a[],double b[],double c[], double f[], double yy[],int n)
{
	int i;
	c [0]= c[0]/b[0];
	yy[0]= f[0]/b[0];
	for( i=1; i<n; i++ ){
		b [i]  =  b[i]- a[i]*c[i-1];
		c [i]  =  c[i]/b[i];
		yy[i]  = (f[i]- a[i]*yy[i-1])/b[i];
	}
	for( i=n-2; i>=0; i-- )
		yy[i]= yy[i] - c[i]*yy[i+1];
	return 1;
}

int GaussPivoting( double **A, double p[], int n, double x[] )
{
	double xi;
	int   maxi, m, i,j, it,ip,jt;
	double temp;
	
	for(i=0;i<n;i++)
	{
		// pivoting process
		maxi=i;
		for( m=i+1; m<n; m++)
		{
			if(fabs(A[m][i])>fabs(A[i][i]))
				maxi=m;
		}
		if( maxi!=i ){
		    // exchange matrix A[maxi], A[i]
			for( m=0; m<n; m++)
			{
				temp=A[maxi][m];
				A[maxi][m]=A[i][m];
				A[i][m]=temp;
			}
			// change right hand side
			temp=p[maxi];
			p[maxi]=p[i];
			p[i]=temp;
		}

		if( fabs(A[i][i])<1.e-20 ){
			printf("matrix is singular");
			return 0;
		}

		// elimination process
		for(it=i+1;it<n;it++)
		{
			xi=A[it][i]/A[i][i];  // factor
			A[it][i]=0;
			p[it] -= xi*p[i];
			for( j=i+1;j<n;j++ )
			{
				A[it][j]-=xi*A[i][j];
			}
		}
	}
	
	// back substitution
    x[n-1]=p[n-1]/A[n-1][n-1];
	for(ip=n-2;ip>=0;ip--)
	{
		double sum=p[ip];
		for(jt=n-1;jt>ip;jt--) 
		{
			sum-=A[ip][jt]*x[jt];
		}
		x[ip]=sum/A[ip][ip];
	}
	return 1;
}


int GaussSeidelIter( double **A, double b[], double x0[], int n, double x[] )
{
//	int i;
	do{
		
	}while(1);
	return 1;
}
