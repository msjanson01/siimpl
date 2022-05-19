/* mtrxfun.c  */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mtrxfun.h"
#include "mjmath.h"

/* 3D vector functions -----------------------------------------------*/
void addeq3D(struct double3D *a, const struct double3D *b)
{
   a->x += b->x;
   a->y += b->y;
   a->z += b->z;
}

struct double3D add3D(const struct double3D *a,
							 const struct double3D *b)
{
struct double3D res;
   res.x = a->x + b->x;
   res.y = a->y + b->y;
   res.z = a->z + b->z;
   return res;
};

void subeq3D(struct double3D *a, const struct double3D *b)
{
   a->x -= b->x;
   a->y -= b->y;
   a->z -= b->z;
}

struct double3D sub3D(const struct double3D *a,
							 const struct double3D *b)
{
struct double3D res;
   res.x = a->x - b->x;
   res.y = a->y - b->y;
   res.z = a->z - b->z;
   return res;
};

double norm3D(const struct double3D *a)
{
  	return sqrt(sqr(a->x) + sqr(a->y) + sqr(a->z));
}

double sum3D(const struct double3D *a)
{
  	return ((a->x) + (a->y) + (a->z));
}

double dot3D(const struct double3D *a,
			  const struct double3D *b)
{
  	return a->x*b->x + a->y*b->y + a->z*b->z;
};

struct double3D cross3D(const struct double3D *a,
				 const struct double3D *b)
{
	struct double3D res;
   res.x = a->y * b->z - a->z * b->y;
   res.y = a->z * b->x - a->x * b->z;
   res.z = a->x * b->y - a->y * b->x;
   return res;
};

struct double3D timesk3D(const struct double3D *a,const double k)
{
	struct double3D res;

   res.x = (a->x)*k;
   res.y = (a->y)*k;
   res.z = (a->z)*k;

   return res;
}

struct double3D divk3D(const struct double3D *a,const double k)
{
	struct double3D res;

   res.x = (a->x)/k;
   res.y = (a->y)/k;
   res.z = (a->z)/k;

   return res;
}


void print3D(const struct double3D a)
{
	printf("%f %f %f\n",a.x,a.y,a.z);
};
void printInt3D(const struct int3D a)
{
	printf("%d %d %d\n",a.x,a.y,a.z);
};

/* -- Functions including base3D structure -----*/

struct double3D vecttransf(struct double3D a1,struct base3D C1)
{
	struct double3D a2;

   a2.x = a1.x*C1.a1.x + a1.y*C1.a2.x + a1.z*C1.a3.x;
   a2.y = a1.x*C1.a1.y + a1.y*C1.a2.y + a1.z*C1.a3.y;
   a2.z = a1.x*C1.a1.z + a1.y*C1.a2.z + a1.z*C1.a3.z;

	return a2;
};

struct double3D floorvecttransf(struct double3D a1,struct base3D C1)
{
	struct double3D a2;

   a2.x = floor(a1.x*C1.a1.x + a1.y*C1.a2.x + a1.z*C1.a3.x);
   a2.y = floor(a1.x*C1.a1.y + a1.y*C1.a2.y + a1.z*C1.a3.y);
   a2.z = floor(a1.x*C1.a1.z + a1.y*C1.a2.z + a1.z*C1.a3.z);

	return a2;
};

struct base3D vects2base(struct double3D a1, struct double3D a2,
								 struct double3D a3)
{
	struct base3D base;
   base.a1 = a1;
   base.a2 = a2;
   base.a3 = a3;
   return base;
};

double volumebase3D(struct base3D B)
{
	struct double3D a;

		a = cross3D( &B.a1, &B.a2);

	return dot3D( &a, &B.a3);
};



/* functions to calculate inverse 3*3 matrix (struct base3D)  */
double subdet(double a,double b,double c,double d)
{
	return a*d-b*c;
};
double detbase(struct base3D A)
{
   return
	A.a1.x*subdet(A.a2.y,A.a2.z,A.a3.y,A.a3.z)-
   A.a1.y*subdet(A.a2.x,A.a2.z,A.a3.x,A.a3.z)+
	A.a1.z*subdet(A.a2.x,A.a2.y,A.a3.x,A.a3.y);
};
struct  base3D invbase(struct base3D A)
{
	struct base3D B;
   double a,b,c,d,e,f,g,h,i,detA;

   a=A.a1.x;b=A.a1.y;c=A.a1.z;
   d=A.a2.x;e=A.a2.y;f=A.a2.z;
   g=A.a3.x;h=A.a3.y;i=A.a3.z;
   detA=detbase(A);

   B.a1.x=subdet(e,f,h,i)/detA;
   B.a2.x=-subdet(d,f,g,i)/detA;
 	B.a3.x=subdet(d,e,g,h)/detA;

   B.a1.y=-subdet(b,c,h,i)/detA;
   B.a2.y=subdet(a,c,g,i)/detA;
 	B.a3.y=-subdet(a,b,g,h)/detA;

   B.a1.z=subdet(b,c,e,f)/detA;
   B.a2.z=-subdet(a,c,d,f)/detA;
 	B.a3.z=subdet(a,b,d,e)/detA;

	return B;
};
/* -----------------------------------------------------------------*/

