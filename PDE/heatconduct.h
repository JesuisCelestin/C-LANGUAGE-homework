#ifndef HEATCONDUCT_H
#define HEATCONDUCT_H

void HeatConduct_Init( int tNX,double ttotaltime, double tromda );
void HeatConduct_FTCS();
void HeatConduct_BTCS();
void HeatConduct_CrankNilson();
void HeatConduct_4th();
void HeatConduct_Output();

#endif
