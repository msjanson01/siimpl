// file: mcstructdef.h
// Martin Janson, 2000-03-10


// mc structures
struct SCstruct {
   int id;
	struct double3D r;
};

struct vibStruct {
	double A;
   struct double3D k;
};

struct ionStruct {
	int id;
	double E;
	struct double3D ev0;
	struct double3D r;
	struct double3D rUC;
	struct double3D UCr0;
   int SCi;
   double w;
   double Rpath;
   struct vibStruct vib;
};

struct atomStruct {
	double x;
   int Z;
   double M;
	double Ed;
	double u3;
};

struct espStruct {
	double A1;
   double p;
   double p2;
   double A3;
   double A4;
   double A5;
   double s;
   double lf;
};

// ill = ion linked list, for recoils and rare event algorithms
struct illStruct {
   struct ionStruct ion;
	struct illStruct *p;
};

// sc = simultaneous collision
struct scStruct {
	int tgatomi;
   struct double3D RcUC;
   struct double3D R2UC;
   double dRv0;
   double p;

};

// -------------------------------------------------------


