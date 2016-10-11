#include <stdio.h>
#include <math.h>
#include "ns2d.h"
#include "array.h"

static const double Re=400.,sorwei=1.8;
static int    nx,ny;
static double **u,**v,**vor,**stf, dx,dy,dx2inv,dy2inv,dsinv,dxinv,dyinv;
static double dt=1.e-3,ttime;

void ns2d_init(int tnx,int tny)
{
	int i,j;
	nx=tnx;
	ny=tny;
	dx= 1./nx;
	dy= 1./ny;
	dx2inv= 1. /(dx*dx);
	dy2inv= 1. /(dy*dy);
	dsinv = 0.5/(dx2inv+dy2inv);
	dxinv= 1./dx;
	dyinv= 1./dy;

	u  = Make2DArray(nx+1,ny+1);
	v  = Make2DArray(nx+1,ny+1);
	vor= Make2DArray(nx+1,ny+1);
	stf= Make2DArray(nx+1,ny+1);

	// init flows
	for(i=0;i<=nx;i++)
	for(j=0;j<=ny;j++)
	{
		u  [i][j]= 0.;
		v  [i][j]= 0.;
		vor[i][j]= 0.;
		stf[i][j]= 0.;
	}
	for(i=1;i<nx;i++)
		u  [i][ny]= 1.;
}

int ns2d_solve()
{
	int n;
	ttime=0.;
	for( n=0;n<500000;n++ )
	{
		ns2d_iter();
		ttime= ttime+dt;
		if( n%50==0 )
			printf("n=%d, time=%f \n", n, ttime);
		if( n%500==0 ){
			ns2d_output();
			// return 0;
		}
	}
	return 1;
}

void ns2d_iter()
{
	int i,j,iter;
	double fxl,fxr,fyl,fyr,flux,fluy,Res,Res1,tmp,ttt;
	// vortex evolution, solve using first order Euler method
	//    boundary
	for(i=1;i<nx;i++) {
		vor[i][0] = dyinv*( u[i][1]  - u[i][0] );
		vor[i][ny]= dyinv*( u[i][ny] - u[i][ny-1] );
	}
	for(j=1;j<ny;j++) {
		vor[0][j] = -dxinv*( v[1] [j] - v[0]   [j] );
		vor[nx][j]= -dxinv*( v[nx][j] - v[nx-1][j] );
	}

	for( i=1;i<nx;i++ )
	for( j=1;j<ny;j++ )
	{
		fxr= ( vor[i+1][j  ] - vor[i  ][j  ] ) * dxinv;
		fxl= ( vor[i  ][j  ] - vor[i-1][j  ] ) * dxinv;
		fyr= ( vor[i  ][j+1] - vor[i  ][j  ] ) * dyinv;
		fyl= ( vor[i  ][j  ] - vor[i  ][j-1] ) * dyinv;
		flux= -0.5*(( u[i][j] - fabs(u[i][j] ))*fxr + (u[i][j] + fabs(u[i][j]))*fxl );
		fluy= -0.5*(( u[i][j] - fabs(u[i][j] ))*fyr + (u[i][j] + fabs(u[i][j]))*fyl );

		vor[i][j]= vor[i][j]+ dt*( flux + fluy +
					1./Re*( dx2inv*( vor[i+1][j]- 2.*vor[i][j]+ vor[i-1][j] ) +
							dy2inv*( vor[i][j+1]- 2.*vor[i][j]+ vor[i][j-1] )  ) );
	}

	// stream function, solving using SOR iteration
	for(iter=0; iter<100; iter++)
	{
		Res = 0.;
		for(i=1;i<nx;i++)
		for(j=1;j<ny;j++)
		{
			tmp = stf[i][j];
			ttt = dsinv*( dx2inv*(stf[i+1][j] + stf[i-1][j]) +
						  dy2inv*(stf[i][j+1] + stf[i][j-1]) - vor[i][j] ) ; // Gaussian-Seidel iteration;
			stf[i][j]= (1.-sorwei)*stf[i][j] + sorwei*ttt ;
			Res += fabs( stf[i][j] - tmp );
		}
		if( iter==0 )
			Res1= Res;
		else if( Res/Res1<1.e-2 )
			break;
	}
	// printf("%d,%f\n",iter,Res/Res1);

	// velocity
	for(i=1;i<nx;i++)
	for(j=1;j<ny;j++)
	{
		u[i][j]= 0.5* dyinv * ( stf[i][j+1] - stf[i][j-1] );
		v[i][j]=-0.5* dxinv * ( stf[i+1][j] - stf[i-1][j] );
	}
}

// output results, tecplot format
void ns2d_output()
{
	int i,j;
	FILE *fp;
	fp= fopen("out_ns.dat","w+");
	fprintf(fp," variables=\"x\",\"y\",\"u\",\"v\",\"vorticity\",\"stream line\" ");
	fprintf(fp," zone i=%d   j=%d f=point ", nx+1,ny+1);
	for(i=0;i<=nx;i++)
	for(j=0;j<=ny;j++)
	{
		fprintf( fp,"%f,%f,%f,%f,%f,%f \n ",dx*i,dy*j,u[i][j],v[i][j],vor[i][j],stf[i][j] );
	}
	fclose(fp);
}
