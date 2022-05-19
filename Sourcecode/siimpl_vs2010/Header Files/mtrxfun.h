/* mtrxfun.h  */


/* Define structures ----------------*/

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


// function prototypes
/* 3D vector functions -----------------------------------------------*/
void addeq3D(struct double3D *a, const struct double3D *b);

struct double3D add3D(const struct double3D *a,
							 const struct double3D *b);
                      
void subeq3D(struct double3D *a, const struct double3D *b);

struct double3D sub3D(const struct double3D *a,
							 const struct double3D *b);

double norm3D(const struct double3D *a);

double sum3D(const struct double3D *a);

double dot3D(const struct double3D *a,
			  	 const struct double3D *b);

struct double3D cross3D(const struct double3D *a,
				 				const struct double3D *b);

struct double3D timesk3D(const struct double3D *a,const double k);

struct double3D divk3D(const struct double3D *a,const double k);

void print3D(const struct double3D a);

void printInt3D(const struct int3D a);


/* -- Functions including base3D structure -----*/

struct double3D vecttransf(struct double3D a1,struct base3D C1);

struct double3D floorvecttransf(struct double3D a1,struct base3D C1);

struct base3D vects2base(struct double3D a1, struct double3D a2,
								 struct double3D a3);

double volumebase3D(struct base3D B);



/* functions to calculate inverse 3*3 matrix (struct base3D)  */
double subdet(double a,double b,double c,double d);

double detbase(struct base3D A);

struct  base3D invbase(struct base3D A);
/* -----------------------------------------------------------------*/


