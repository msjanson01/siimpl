// file: mcglobalextern.h
// Martin Janson, 2000-03-09

// --- Globaly used definitions
#define MAXSCATOMS 100
#define MAXNSIMCOLL 10

// BCA
// DRMINMIN is the smallest advance for an ion not to
// be considdered as a self collision. i.e. EPS
#define DRMINMIN 1e-10

#define MAXNRAREDEPTH 5

//  IAPfunction codes
#define MAGICFORMULA 0
#define GAUSSMEHLERZBL 1

// Thermal vibration algorithm codes
#define ROBINSONVIB 0
#define MJVIB 1

// Constants
#define INF 1e100

// extern GLOBAL varaibles defined in mcglobaldef -----------------------

// Crystal structure data  -----------
// target data
extern int Nspecies;
extern struct atomStruct *atom;
extern double nAtom, nAtomAmorph;
extern int randomtarget;

// Unit Cell
extern struct base3D UCbase,on2UCbase;
extern struct SCstruct *UCatom;
extern int NatUC; // Number of atoms in unit cell


// Sub Cells
extern struct base3D SCbase,on2SCbase;
extern struct int3D nSC;
extern int NatSC, NSC;
extern double SCsize;
extern struct SCstruct *SCatom;
extern int 	 *SCpointer;

// ----------------------------------

// Primary implanted ions data ---------
extern struct ionStruct ion0;
extern double dose;

// ----------------------------------

// target surface and layer width  -----
extern struct double3D Rsurf;
extern struct double3D ensurf, e1surf, e2surf;
extern double layerWidth;
extern double randsurflayer;
extern double U_s; // (keV)


//  implanted area definition ----------
extern struct double3D Rimp;
extern double impArea1, impArea2;

// Nuclear stopping data
extern int IAPfunction;

// electronic stopping power data
extern struct espStruct *esp;
extern double E_fermi;
extern double *locnorm;

// BCA paramters
extern double pmax,ppmax,Estop;
extern double dRmin, dRsc; // (Å)
extern int simcollalg, vibalg;

// Rare event variables
extern int rareEventMult;
extern double rareDepth[MAXNRAREDEPTH + 1];
extern int NrareDepth;

//	-----   ion, damage, energy, distribution variables -------
extern struct illStruct *ill;


// Damage model
extern int recoil;
extern double rndscatfac;

// input statistics variables ---------------------
extern int	Nions;
extern int Ndistz;
extern double dz;

// output file identifyers (should later not be defined globally!!)
extern FILE *fionR3, *fIR3, *fVR3;

extern double *Czion, *CzI, *CzV;
extern int saveR3;
extern int saveionsR3;
extern int saverecoilsR3;

extern double savedose[10];
extern int Nsavedose;


extern double *Qa, *Ta;

extern int calcFnep;
extern double *Fnep;

// miscellaneous
extern long randomz;

// -------------------------------------------------------

