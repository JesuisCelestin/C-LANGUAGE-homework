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
int m=5;
int k=5;
int n=5;
int h=3;


void HeatConduct_Init( int tNX,double ttotaltime, double tromda )   //初始化
{
    int i;
    double xx;
    romda    = tromda;
    totaltime= ttotaltime;
    NX       = tNX;
    //
    dx = 1./NX;    //对x划分
    dt = 4.e-5;   // 0.8 * dx*dx*0.5/romda;   // change this to given physical time step in implicit scheme
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
		u[0] = 0.;    ///u(0,t)=0, n(1,t)=0;
		u[NX]= 0.;
		// inner points
		for(i=1; i<NX; i++){
			u[i]= sigma*(u1[i-1] + u1[i+1]) + (1.-2*sigma)*u1[i];
		}    ///u[i]储存Tn+1[i] u1储存Tn[i]
		// store backup
		for(i=0; i<=NX; i++)
			u1[i]= u[i];   ///刷新

		if(istop==1)break;
	}
     ///因为示例代码的output函数使得运行程序一次只能输出一个heat.dat文件，不方便画图，故作出以下修改
	//打印输出
	char str[]="heat_FTCS0.X.dat";
	str[11]=m+'0';

	FILE *fp;
	fp = fopen(str,"w+");
	for(i=0;i<=NX;i++)
		fprintf(fp,"%f,%f\n",i*dx,u[i]);
	fclose(fp);
    m--;

}

void HeatConduct_BTCS()
{
	int istep,i,j,MaxCycle, istop=0;
	double sigma= romda*dt/(dx*dx), timep, *av,*bv,*cv,*fv;

	MaxCycle= 1000000;
	timep   = 0.;
	printf("sigma:%f\n",sigma);
	av  = Make1DArray(NX-1);
	bv  = Make1DArray(NX-1);
	cv  = Make1DArray(NX-1);
	fv  = Make1DArray(NX-1);



    for(istep=0; istep<MaxCycle; istep++){
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = romda*dt/(dx*dx);
			istop = 1;
		}
		// bounary condition
		u[0] = 0.;    ///u(0,t)=0, n(1,t)=0;
		u[NX]= 0.;
		// inner points
	    for(i=0;i<NX-1;i++)
        {
            bv[i]=1+2*sigma;
            av[i]=-sigma;
            cv[i]=-sigma;
            fv[i]=u1[i+1];
        }

        chase(av,bv,cv,fv,u+1,NX-1);


		///对一个n, 解出了所有的T[i]
		// store backup

       for(i=1;i<=NX-1;i++)
         {
             u1[i]=u[i];
         }

		   if(istop==1)break;
	}
	free(av);
	free(bv);
	free(cv);
	free(fv);


    char str[]="heat_BTCS0.X.dat";
	str[11]=n+'0';

    FILE *fp;
    fp = fopen(str,"w+");
	for(i=0;i<=NX;i++)
		fprintf(fp,"%f,%f\n",i*dx,u[i]);
	fclose(fp);
	n--;

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

   for(istep=0; istep<MaxCycle; istep++){
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = romda*dt/(2*dx*dx);
			istop = 1;
		}
		// bounary condition
		u[0] = 0.;    ///u(0,t)=0, n(1,t)=0;
		u[NX]= 0.;
		// inner points
	    for(i=0;i<NX-1;i++)
        {
            bv[i]=2+(1/sigma);
            av[i]=-1;
            cv[i]=-1;
            fv[i]=u1[i+2]-(2-1/sigma)*u1[i+1]+u1[i];
        }

        chase(av,bv,cv,fv,u+1,NX-1);


		///对一个n, 解出了所有的T[i]
		// store backup

       for(i=1;i<=NX-1;i++)
         {
             u1[i]=u[i];
         }



		   if(istop==1)break;
	}

	free(av);
	free(bv);
	free(cv);
	free(fv);

	char str[]="heat_CrankNilson.X.dat";
	str[17]=k+'0';

	FILE *fp;
	fp = fopen(str,"w+");
	for(i=0;i<=NX;i++)
		fprintf(fp,"%f,%f\n",i*dx,u[i]);
	fclose(fp);
    k--;
}


void HeatConduct_4th()
{
	int istep,i,MaxCycle,j,istop=0,iter;
	double sigma= 24.*dx*dx/(dt*romda),sigma1=romda*dt/(dx*dx), timep, *f, errmax,tu,
		wei=1.8;
	MaxCycle= 1000000;
	timep   = 0.;
	printf("sigma:%f\n",sigma);

	f = Make1DArray(NX-1);
	 for(istep=0; istep<MaxCycle; istep++){
		timep += dt;
		if(timep>=totaltime){
			dt= totaltime - (timep-dt);
			timep = totaltime;
			sigma = 24.*dx*dx/(dt*romda);
			istop = 1;
		}


     u[0]=0;
     f[0]=u[1]= sigma1*(u1[2])+(1-2*sigma1)*u1[1]; //u[0]=0

     for(i=1;i<NX-2;i++)
        {
            f[i]=-u1[i+3]+16*u1[i+2]-(30-sigma)*u1[i+1]+16*u1[i]-u1[i-1];
        }
        f[NX-2]=u1[NX-1];

    for(iter=0; iter<MaxCycle; iter++)
		{
			errmax= 0.;
			for(i=2; i<NX-1; i++){
				tu = u[i];
				u[i]= (f[i-1]-(u[i-2]-16.*u[i-1]-16.*u1[i+1]+u1[i+2]))/(30+sigma);
				u[i]= (1-wei)*tu+ wei*u[i];
				errmax= max(errmax, fabs(tu-u[i]));
			}
			if( errmax<1.e-5 ) break;
		}
        u[NX-1]= sigma1*(u1[NX-2])+(1-2*sigma1)*u1[NX-1];
        u[NX]=0;

		for(i=0;i<=NX;i++)
            u1[i]=u[i];

      if(istop==1)break;
	 }

	free(f);
	FILE *fp;

	fp = fopen("HeatConduct_4th_with_NX40.dat","w+");
	for(i=0;i<=NX;i++)
		fprintf(fp,"%f,%f\n",i*dx,u[i]);
	fclose(fp);
}


void HeatConduct_Output()
{
	int i;
	FILE *fp;
	fp = fopen("heat.dat","w+");
	for(i=0;i<=NX;i++)
		fprintf(fp,"%f,%f\n",i*dx,u[i]);
	fclose(fp);

}
