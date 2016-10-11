#include <math.h>
#include <stdio.h>
#include "heatconduct2d.h"
#include "linear.h"
#include "array.h"

#define pi        3.14159265358
#define max(x,y)  x>y?x:y

// global variables
static int NX,NY;
static double **u, **u1, dx,dy,dt, totaltime,romda;

void HeatConduct2D_Init( int tNX,int tNY, double ttotaltime, double tromda )
{
	int i,j;
    double xx,yy;
    romda    = tromda;
    totaltime= ttotaltime;
    NX       = tNX;
	NY       = tNY;
    //
    dx = 1./NX;
	dy = 1./NY;
    dt = 0.001;   // 0.8 * dx*dx*0.5/romda;   // change this to given physical time step in implicit scheme
	u1 = Make2DArray(NX+1,NY+1);
	u  = Make2DArray(NX+1,NY+1);
    for(i=0;i<=NX;i++)
	for(j=0;j<=NY;j++){
    	xx   = dx*i;
		yy   = dy*j;
        u [i][j] = exp( -20.*(xx-0.5)*(xx-0.5) -20.*(yy-0.5)*(yy-0.5) );
        u1[i][j] = u[i][j];
    }
}

void HeatConduct2D_SOR()
{
	int i,j,istep,istop=0,iter;
	double timep,sigmax,sigmay,sigmainv, **f,tu,sorwei=1.0,errmax;

	sigmax = romda*dt/(2.*dx*dx);
	sigmay = romda*dt/(2.*dy*dy);
	sigmainv = 1./(1.+2*sigmax+2*sigmay);
	f= Make2DArray(NX+1,NY+1);
	timep=0.;

	// boundary
	for(i=0;i<=NX;i++){
		u[i][0 ]= 0.;
		u[i][NY]= 0.;
	}
	for(j=0;j<=NY;j++){
		u[0 ][j]=0.;
		u[NX][j]=0.;
	}

	for(istep=0; istep<10000000; istep++)
	{
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigmax = romda*dt/(2.*dx*dx);
			sigmay = romda*dt/(2.*dy*dy);
			sigmainv = 1./(1.+2*sigmax+2*sigmay);
			istop = 1;
		}
		printf("time: %f \n",timep);
		// right hand side
		for(i=1; i<NX; i++)
		for(j=1; j<NY; j++){
				f[i][j]= sigmax*(u[i+1][j]-2.*u[i][j]+u[i-1][j]) +
					     sigmay*(u[i][j+1]-2.*u[i][j]+u[i][j-1]) + u[i][j];
		}
		// inner points, SOR iteraction
		for(iter=0; iter<100; iter++)
		{
			errmax= 0.;
			for(i=1; i<NX; i++)
			for(j=1; j<NY; j++){
				tu = u[i][j];
				u[i][j]= sigmainv*( f[i][j] +sigmax*(u[i+1][j]+u[i-1][j]) +
					                        +sigmay*(u[i][j+1]+u[i][j-1]) );
				u[i][j]= (1.-sorwei)*tu+ sorwei*u[i][j];
				errmax= max( errmax, fabs(tu-u[i][j]) );
			}
			if( errmax<1.e-5 ) break;
		}
		//HeatConduct2D_Output();return;
		printf("iter=%d, errmax=%f\n",iter,errmax);

		if(istop==1)break;
	}

	free2DArray(f,NX+1);
}

void HeatConduct2D_ADI()
{
}


void HeatConduct2D_Output()
{
	int i,j;
	double tu;
	FILE *fp;
	fp= fopen("heat2d.dat","w+");
	fprintf(fp," variables=\"x\",\"y\",\"u\"\n ");
	fprintf(fp," zone i=%d   j=%d f=point \n ", NX+1,NY+1);
	for(i=0;i<=NX;i++){
        for(j=0;j<=NY;j++){
		tu = u[i][j];
		fprintf(fp,"%f  ",tu);
	       }
        fprintf(fp,"\n");
	}

	fclose(fp);
	free2DArray(u, NX+1);
	free2DArray(u1,NX+1);
}
