#include <math.h>
#include <stdio.h>
#include "heatconduct.h"
#include "linear.h"
#include "array.h"

#define pi        3.14159265358
#define max(x,y)  x>y?x:y

// global variables
static int NX;
static double *u, *u1, dx,dt, totaltime,romda;

void HeatConduct_Init( int tNX,double ttotaltime, double tromda )
{
    int i;
    double xx;
    romda    = tromda;
    totaltime= ttotaltime;
    NX       = tNX;
    //
    dx = 1./NX;
    dt = 0.00001;   // 0.8 * dx*dx*0.5/romda;   // change this to given physical time step in implicit scheme
    u  = Make1DArray(NX+1);
    u1 = Make1DArray(NX+1);
    for(i=0;i<=NX;i++){
    	xx   = dx*i;
        u [i] = sin(pi*xx) + sin(3*pi*xx);
        u1[i] = u[i];
    }
}

void HeatConduct_FTCS()
{
	int istep,i,MaxCycle, istop= 0;
	double sigma= romda*dt/(dx*dx), timep;
	MaxCycle= 1000000;
	timep   = 0.;
	printf("sigma:%f\n",sigma);
	for(istep=0; istep<MaxCycle; istep++){
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = romda*dt/(dx*dx);
			istop = 1;
		}
		// bounary condition
		u[0] = 0.;
		u[NX]= 0.;
		// inner points
		for(i=1; i<NX; i++){
			u[i]= sigma*(u1[i-1] + u1[i+1]) + (1.-2*sigma)*u1[i];
		}
		// store backup
		for(i=0; i<=NX; i++)
			u1[i]= u[i];

		if(istop==1)break;
	}
}

void HeatConduct_BTCS()
{
	int istep,i,MaxCycle,j, istop=0;
	double sigma= romda*dt/(dx*dx), timep, *av,*bv,*cv,*fv;
	MaxCycle= 1000000;
	timep   = 0.;
	printf("sigma:%f\n",sigma);
	av  = Make1DArray(NX-1);
	bv  = Make1DArray(NX-1);
	cv  = Make1DArray(NX-1);
	fv  = Make1DArray(NX-1);


for(istep=0; istep<=MaxCycle; istep++){
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = romda*dt/(dx*dx);
			istop = 1;
		}
		//边界条件
		u[0]=0;
		u[NX]=0;
		for (i=0;i<(NX-1);i++)
    { av[i]=-sigma;
      bv[i]=2*sigma+1;
      cv[i]=-sigma;
      fv[i]=u[i+1];
  }                          //定义chase法矩阵数组


		chase(av,bv,cv,fv,&u[1],NX-1);


	if(istop==1)break;
    }
//for(i=0;i<(NX-1);i++)
    //{u[i+1]=yy[i];}
	free(av);
	free(bv);
	free(cv);
	free(fv);
}

void HeatConduct_CrankNilson()
{
	int istep,i,MaxCycle,j,istop=0;
	double sigma= romda*dt/(2.*dx*dx), timep, *av,*bv,*cv,*fv;
	MaxCycle= 1000000;
	timep   = 0.;
	printf("sigma:%f\n",sigma);
	av  = Make1DArray(NX-1);
	bv  = Make1DArray(NX-1);
	cv  = Make1DArray(NX-1);
	fv  = Make1DArray(NX-1);

for(istep=0; istep<=MaxCycle; istep++){
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = romda*dt/(dx*dx);
			istop = 1;
		}

    u[0]=0;
   u[NX]=0;
		for (i=0;i<(NX-1);i++)
    { av[i]=-1;
      bv[i]=(2+(1/sigma));
      cv[i]=-1;
      fv[i]=-((2-(1/sigma))*u[i+1])+u[i]+u[i+2];
  }                          //定义chase法矩阵数组


		chase(av,bv,cv,fv,&u[1],NX-1);


	if(istop==1)break;
    }

	free(av);
	free(bv);
	free(cv);
	free(fv);
}

void HeatConduct_4th()
{
	int istep,i,MaxCycle,j,istop=0,iter;
	double sigma0= romda*dt/(2*dx*dx),sigma=(24*dx*dx)/(romda*dt),timep, *b, errmax,tu,
		wei=1.8;
	MaxCycle= 1000000;
	timep   = 0.;
	iter    = 0;
	printf("sigma:%f\n sigma0:%f\n",sigma,sigma0);
	b = Make1DArray(NX-1);




	for(istep=0; istep<=MaxCycle; istep++){

		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = romda*dt/(dx*dx);
			istop = 1;
		}

     //给出初值。
      u[0]=0;
      u[NX]=0;




     //定义矩阵右端项
        b[0]=(-(2-(1/sigma0)))*u[1]+u[2];
                                               //u[1]、u[NX-1]使用二阶显式格式
     for (i=1;i<(NX-2);i++)
     {
         b[i]=-u[i+3]+16*u[i+2]-(30-sigma)*u[i+1]+16*u[i]-u[i-1];
     }                                             //u[2]~u[NX-2]使用四阶格式
      b[NX-2]=(-(2-(1/sigma0)))*u[NX-1]+u[NX-2];


     //u的迭代

     do{
         errmax=0.;
	         tu=u[1];
	       u[1]=1/(2+1/sigma0)*(b[0]+u[2]);
	       u[1]=(1-wei)*tu+wei*u[1];
         errmax= max( errmax, fabs(tu-u[1]) );//u[1]单独拿出来


	      for(i=2;i<NX-1;i++)
          {
              tu=u[i];
              u[i]=1/(30+sigma)*(b[i-1]+16*u[i-1]-u[i-2]+16*u[i+1]-u[i+2]);
              u[i]=(1-wei)*tu+wei*u[i];

              errmax= max( errmax, fabs(tu-u[i]) );
          }

               tu=u[NX-1];
	      u[NX-1]=1/(2+1/sigma0)*(b[NX-2]+u[NX-2]);
	      u[NX-1]=(1-wei)*tu+wei*u[NX-1];
	       errmax= max( errmax, fabs(tu-u[NX-1]) );//u[NX-1]单独拿出来



        }while(errmax>1.e-6);


       // printf("iter=%d, errmax=%f\n",iter,errmax);










	if(istop==1)break;
}


	free(b);

}

void HeatConduct_Output()
{
	int i;
	FILE *fp;
	fp = fopen("heat.dat","w+");
	for(i=0;i<=NX;i++)
		fprintf(fp,"%f,%f\n",i*dx,fabs(u[i]));
	fclose(fp);
	printf("lalalala");
}

