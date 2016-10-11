#include <stdio.h>
#include <math.h>
#include "array.h"
#include "poisson.h"

#define pi 3.14159265359

static const double sorwei=1.9;
static int    nx, ny;
static double **u, dx,dy;

void poisson_init( int tnx, int tny )
{
	int i,j;
	double xx,yy;

	nx= tnx;
	ny= tny;
	dx= 1./nx;
	dy= 1./ny;
	u= Make2DArray(nx+1,ny+1);

	// init field
	for(i=0;i<=nx;i++)
	for(j=0;j<=ny;j++)
	{
		xx = i*dx;
		yy = j*dy;
		u[i][j] = 1.; // 1+xx+sin(pi*xx)*sin(pi*yy);
	}
}

int poisson_sor( )
{
	int i,j, iter=0;
	double dx2inv,dy2inv,dsinv, Res,Res1,ttt,tu,sss;

	dx2inv= 1. /(dx*dx);
	dy2inv= 1. /(dy*dy);
	dsinv = 0.5/(dx2inv+dy2inv);



	return 1;
}

int poisson_4th_sor( )
{
	int i,j, iter=0;
	double dx2inv,dy2inv,dsinv, dx12,dy12,ds25, Res,Res1,ttt,tu,sss;

	dx2inv= 1. /(dx*dx);
	dy2inv= 1. /(dy*dy);
	dsinv = 0.5/(dx2inv+dy2inv);

	dx12 = 1. /12./(dx*dx);
	dy12 = 1. /12./(dy*dy);
	ds25 = 0.4/(dx2inv+dy2inv);

	// form five-diangonal matrix
	for(iter=0; iter<500000; iter++)
	{
		// boundary
		for( i=0;i<=nx;i++ )
		{
			u[i][0] = 1.+i*dx;
			u[i][ny]= 1.+i*dx;
		}
		for( j=0;j<=ny;j++ )
		{
			u[0 ][j]= 1.;
			u[nx][j]= 1.+nx*dx;
		}

		// inner points
		Res = 0.;
		for( i=1; i<nx; i+=nx-2)
		for( j=1; j<ny; j++){
			ttt= u[i][j];
			sss= -2*pi*pi*sin(pi*i*dx)*sin(pi*j*dy);  // source term
			tu = dsinv*( dx2inv*(u[i+1][j]+u[i-1][j]) +
						 dy2inv*(u[i][j+1]+u[i][j-1]) - sss ) ; // Gaussian-Seidel iteration
			u[i][j]= (1.-sorwei)*u[i][j] + sorwei*tu;
			Res += fabs( u[i][j]-ttt );
		}
		for( i=1; i<nx; i++     )
		for( j=1; j<ny; j+=ny-2 ){
			ttt= u[i][j];
			sss= -2*pi*pi*sin(pi*i*dx)*sin(pi*j*dy);  // source term
			tu = dsinv*( dx2inv*(u[i+1][j]+u[i-1][j]) +
						 dy2inv*(u[i][j+1]+u[i][j-1]) - sss ) ; // Gaussian-Seidel iteration
			u[i][j]= (1.-sorwei)*u[i][j] + sorwei*tu;
			Res += fabs( u[i][j]-ttt );
		}


		for( i=2; i<nx-1; i++ )
		for( j=2; j<ny-1; j++ ){
			ttt= u[i][j];
			sss= -2*pi*pi*sin(pi*i*dx)*sin(pi*j*dy);  // source term
			tu =  ds25*( dx12*( -u[i+2][j]-u[i-2][j]+16.*(u[i+1][j]+u[i-1][j]) ) +
				         dy12*( -u[i][j+2]-u[i][j-2]+16.*(u[i][j+1]+u[i][j-1]) ) - sss );
				// Gaussian-Seidel iteration
			u[i][j]= (1.-sorwei)*u[i][j] + sorwei*tu;
			Res += fabs( u[i][j]-ttt );
		}

		if( iter==0 )
			Res1= Res;
		else if( Res/Res1<1.e-12 )
			break;
		if( iter%1000==0 )
		{
			printf("iter=%d, Res=%f \n",iter,Res);
			poisson_output();
		}
	}

	return 1;
}

void poisson_output()
{
	int i,j;
	double x,y,errL1;
	FILE *fp;
	fp= fopen("out_poisson.dat","w+");
	fprintf(fp," variables=\"x\",\"y\",\"u\" \n");
	fprintf(fp," zone i=%d   j=%d f=point ", nx+1,ny+1);
	for(i=0;i<=nx;i++)
	for(j=0;j<=ny;j++)
	{
		fprintf(fp,"%f,%f,%f \n",dx*i,dy*j,u[i][j] );
	}
	fclose(fp);

	errL1 = 0.;
	for( i=1;i<nx;i++ )
	for( j=1;j<ny;j++ )
	{
		x= i*dx; y=j*dy;
		errL1 += fabs( u[i][j] - (1.+x+sin(pi*x)*sin(pi*y)) );
	}
	printf("errL1=%16.12f \n ", errL1);
}
