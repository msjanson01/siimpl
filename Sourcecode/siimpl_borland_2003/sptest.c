// spstd.c

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "spstd.h"

#include "mjmath.h"
#include "mtrxfun.h"
#include "physconstCGS.h"

void rand3D(struct double3D *ev)
// creates random 3D direction (courtecy of UZ)
{
	double phi, Theta;

	phi = ( (double)rand() / RAND_MAX)*2*pi;
	Theta = acos(2*( (double)rand() / RAND_MAX )-1);

   printf("phi = %e, Theta = %e\n", phi, Theta);
   
   ev->x = sin(Theta)*cos(phi);
   ev->y = sin(Theta)*sin(phi);
   ev->z = cos(Theta);
}

void newrandomseed(void)
{
	randomize();
}

int main(void)
{
   double M1,M2,E,p;
   int Z1,Z2;
   int i;

   struct double3D r3;

   Z1 = 10 ; M1 = 20 ;
   Z2 = 1 ; M2 = 2 ;
   E = 50;
   p = 1;

   newrandomseed();
   printf("Z1 = %d, M1 = %f, Z2 = %d, M2 = %f\n",Z1,M1,Z2,M2);
   printf("E = %f, p = %f\n\n",E,p);

	printf("Theta magicFormula = %e\n",
   	2*asin(sqrt(magicFormula( keV2epsl(E,Z1,M1,Z2,M2), 1.0/aZBL(Z1,Z2) ))) );
	printf("ThetaGM = %e\n\n",
   	ThetaGM( keV2epsl(E,Z1,M1,Z2,M2), 1.0/aZBL(Z1,Z2) ) );

	printf("SeLS = %e\n",SeLS(E,Z1,M1,Z2));

   for (i=1; i<=2; i++) {
    rand3D(&r3);
    printf("r3 = %e  %e  %e\n",r3.x,r3.y,r3.z);
   }
   
return 1;
}











