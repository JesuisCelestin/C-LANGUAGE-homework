#ifndef HEATCONDUCT2D_H
#define HEATCONDUCT2D_H

void HeatConduct2D_Init( int tNX,int tNY,double ttotaltime, double tromda );
void HeatConduct2D_SOR();
void HeatConduct2D_ADI();
void HeatConduct2D_Output();

#endif
