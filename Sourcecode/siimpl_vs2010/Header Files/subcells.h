// filename: subcells.h

// function prototypes

void cr_subcells( void);

int insideSphere(struct double3D Rp, struct double3D R0,
	double radius);

int insideCyl(struct double3D Rp, struct double3D R0,
	struct double3D Rcent, double radius);

int insideParep(struct double3D Rp, struct double3D R0,
	struct base3D a);

   