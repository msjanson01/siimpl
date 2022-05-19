// impsim.h


// Define physical constants
#define pi 3.14159265

// Define array size constants
#define MAXSPECIESINUC 2
#define MAXSCATOMS 1100
#define MAXNSC 25



// imsim structures
struct SCstruct {
	struct double3D r;
   int id;
};
struct ionStruct {
   int Z;
   int M;
	double E;
	struct double3D ev0;
	struct double3D r;
	struct double3D rUC;
	struct double3D UCr0;
   int SCi;
   double w;
};
struct atomStruct {
   int Z;
   int M;
	double Ed;
	double u1;
};

// Define GLOBAL varaibles -----------------------

// Crystal structure data  -----------
// Unit Cell
struct base3D UCbase,on2UCbase;

// Sub Cells
struct base3D SCbase,on2SCbase;
struct int3D nSC;
struct SCstruct SCatom[MAXSCATOMS];
int 	 SCpointer[MAXNSC];
// ----------------------------------

// ion and target atom data ---------
struct ionStruct ion0;
struct atomStruct atom[MAXSPECIESINUC];

// ----------------------------------

