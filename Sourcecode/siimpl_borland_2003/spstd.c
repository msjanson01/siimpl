// File: spstd.c
// Martin Janson, 2000-03-21

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "spstd.h"

#include "mjmath.h"
#include "physconstCGS.h"


// Nuclear Stopping Power functions -------------------------------------


// Dimension transfer functions

double keV2epsl(double E0, int Z1, double M1, int Z2, double M2)
{
	return 32.5367*M2*E0/Z1/Z2/(M1+M2)/(pow(Z1,.23)+pow(Z2,.23));
   // 32.5367 = 0.8854*a0(cm)/q^2(esu)*1eV(erg)*1000
}

/* double epsl2keV(double epsl, int Z1, double M1,int Z2,double M2)
{
	return (M1+M2)/M2*Z1*Z2*(pow(Z1,.23)+pow(Z2,.23))/32.5367*epsl;
   // 32.5367 = 0.8854*a0(cm)/q^2(esu)*1eV(erg)*1000
}  */


// Screening constants a
double aZBL(int Z1, int Z2)
{
	return 4.6853e-1/(pow(Z1,.23)+pow(Z2,.23));
	// 4.6853e-1 = 0.8854*a0 (Å)
}
/* double aLin(int Z1, int Z2)
{
	return 4.6853e-1/sqrt(pow(Z1,0.6667)+pow(Z2,0.6667));
	// 4.6853e-1 = 0.8854*a0 (Å)
} */

// ZBL nsp functions
double FiZBL(double x)
{
	return .18175*exp(-3.1998*x)+.50986*exp(-.94229*x)
   			+.28022*exp(-0.4029*x)+.028171*exp(-.20162*x);
}
double grdFiZBL(double x)
{
	return -0.581564*exp(-3.1998*x)-0.480436*exp(-.94229*x)
			 -0.112901*exp(-0.4029*x)-0.00567984*exp(-.20162*x);
}


// find apsis of collision with ZBL potential(dimensionless units)
double apsisx0(double epsl,double b)
{                             
	#define ERRMIN 1e-3
   #define X0START 5
   double x0, x02;

	x0 = X0START;

	do
   {
   	x02 = x0;
     	x0 = -( x0 - sqr(b)/x0 -(.18175*exp(-3.1998*x0)+.50986*exp(-.94229*x0)
                              +.28022*exp(-0.4029*x0)+.028171*exp(-.20162*x0)
                              )/epsl
   			 ) / (1.0 - (-0.581564*exp(-3.1998*x0)-0.480436*exp(-.94229*x0)
			              -0.112901*exp(-0.4029*x0)-0.00567984*exp(-.20162*x0)
                      )/epsl
                  + sqr(b/x0)
                 )
             + x0;

	}  while (fabs(x02/x0 - 1.0) > ERRMIN );

	return x0;

} // end of: apsisx0       -----------------------------



// magicFormula returnes sin^2(Theta/2)
double magicFormula(double epsl, double b)
{

	#define C1 0.99229
	#define C2 0.011615
	#define C3 0.007122
	#define C4 14.813
	#define C5 9.3066

   double sqrtepsl, A, x0, Fi, Rcurv, ans;

	sqrtepsl=sqrt(epsl);

	A = 2.0*(1.0+C1/sqrtepsl)*epsl*pow(b,((C2+sqrtepsl)/(C3+sqrtepsl)));

	x0=apsisx0(epsl,b);

	Fi = FiZBL(x0);

	Rcurv = 2.0*(Fi-epsl*x0)/(grdFiZBL(x0)-Fi/x0);

   // cos(Theta/2)^2 = (b+Rcurv+Delta)/(x0+Rcurv);
	// ans = 1 - cos(Theta/2)^2 = sin(Theta/2)^2;
   ans =  1.0-sqr( (b + Rcurv +
   				A*(x0-b)/(1.0+(C4+epsl)/(C5+epsl)/(sqrt(1.0+sqr(A))-A))
      			)/(x0 + Rcurv) ) ;

   if (ans<=0.0) return 0;
   else return ans;
   
}  // end of: magicFormaula --------------------------



// solution to scatter intergral using Gauss-Mehler quadrature
double ThetaGMzbl(double epsl, double b)
{
	#define NITER 100
	double x0,u;
 	double Hsum = 0.0;
	int j,N_2;


	x0=apsisx0(epsl,b);

	N_2=NITER/2.0;
   for (j=1; j<=N_2; j++) {

		u=cos( pi* ((2.0*j-1.0))/2.0/NITER );
		Hsum += sqrt((1.0-sqr(u))/(1.0-FiZBL(x0/u)*u/x0/epsl-sqr(b*u/x0)));
   }

	return pi-2.0*pi*b/(NITER*x0)*Hsum;

} // end of: ThetaGM ---------------------------





// Electronic Stopping Power Functions  ----------------------------

double SeLS(double E, int Z1, double M1, int Z2)
// Lindhard Sharff Stopping cross section in [keV*Å^2]
{
	// v = sqrt(2*(E*1e3*eV)/(M1*amu)) (cm/s)
   //   = k * sqrt(E(keV)/M1(amu)) (cm/s)
   // 3.8454e-18 = k * (8*pi*eq^2*a0/v0/eV/1e3) [keV*cm^2]

	return 3.8454e-2 * sqrt(E/M1) * pow(Z1,0.166667)*Z1*Z2
   				/pow((pow(Z1,0.666667)+pow(Z2,0.666667)), 1.5);
   // answer in (keV*Å^2)
}











