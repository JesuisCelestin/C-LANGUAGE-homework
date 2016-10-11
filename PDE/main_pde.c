#include <stdio.h>
#include "poisson.h"
#include "ns2d.h"
#include "heatconduct.h"
#include "heatconduct2D.h"

int main( int argc, char *argv[] )
{
	// poisson equations
/*	poisson_init(100,100);
	poisson_4th_sor();
	poisson_output();
*/

	// 2d incompressible navier-stokes equations
	/* ns2d_init(300,100);
	ns2d_solve(); */



	// heat conduction 1D
	/*printf("Using FTCE method\n");
	HeatConduct_Init(100,0.5,1.);
    HeatConduct_FTCS();
    HeatConduct_Init(100,0.4,1.);
    HeatConduct_FTCS();
    HeatConduct_Init(100,0.3,1.);
    HeatConduct_FTCS();
    HeatConduct_Init(100,0.2,1.);
    HeatConduct_FTCS();
    HeatConduct_Init(100,0.1,1.);
    HeatConduct_FTCS();
	//HeatConduct_Output();
	printf("Using BTCE method\n");
    HeatConduct_Init(100,0.5,1.);
	HeatConduct_BTCS();
    HeatConduct_Init(100,0.4,1.);
    HeatConduct_BTCS();
    HeatConduct_Init(100,0.3,1.);
    HeatConduct_BTCS();
    HeatConduct_Init(100,0.2,1.);
    HeatConduct_BTCS();
    HeatConduct_Init(100,0.1,1.);
    HeatConduct_BTCS();
    printf("Using CrankNilson method\n");
    HeatConduct_Init(100,0.5,1.);
    HeatConduct_CrankNilson();
    HeatConduct_Init(100,0.4,1.);
    HeatConduct_CrankNilson();
    HeatConduct_Init(100,0.3,1.);
    HeatConduct_CrankNilson();
    HeatConduct_Init(100,0.2,1.);
    HeatConduct_CrankNilson();
    HeatConduct_Init(100,0.1,1.);
    HeatConduct_CrankNilson();
    HeatConduct_Init(40,0.5,1.);
	HeatConduct_4th();
    HeatConduct_Init(40,0.5,1.);
    HeatConduct_4th();
    HeatConduct_Output(); */





	// heat conduction 2D
	HeatConduct2D_Init(1000,1000,0.1,1.0);
	HeatConduct2D_SOR();
	HeatConduct2D_Output();

	return 0;
}
