/* impsim.c  */

#include <math.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"

/* Define physical constants (CGS) */
#define pi 3.14159265

/* Define array size constants */
#define maxSCatoms 1100
#define maxnSC 25
#define maxFileNameLength 50


/* Define structures */

/* 3D structures */
struct double3D {
	double x;
	double y;
	double z;
};

struct int3D {
	int x;
	int y;
	int z;
};

struct base3D {
	struct double3D a1;
	struct double3D a2;
	struct double3D a3;
};


/* float math. functions --------------------------------------------*/

double sqr(double x)
{
	return x*x;
};

int floormj(double x)
{
	 if (x>=0) return floor(x);
    else return -ceil(-x);
};


void print3D(struct double3D a)
{
	mexPrintf("%f %f %f\n",a.x,a.y,a.z);
};







/* MATLAB interface function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* test variables */
double a,b,c,d,e,f,*in;
int i,j,k,l,m,n;
/* ----------------- '/(

   /* Assign pointers to in- and output arguments */
    in=mxGetPr(prhs[0]);

    a=in[0];
   /* -----------------------------------------------------------------------*/


/* TEST TEST TEST TEST TEST TEST   ------------------------------*/
/*   a.x=1;
   a.y=2;
   a.z=3;

   b.x=4;
   b.y=5;
   b.z=6;
 */
   i=1;
   j=-6;
   a=i/j;
   b=(double) i/j;
   c=(double) (i/j);

   mexPrintf("i/j = %e\n",a);
   mexPrintf("(double) i/j = %e\n",b);
   mexPrintf("(double) (i/j) = %e\n",c);


   return 0;
}
