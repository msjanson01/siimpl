// file: mcglobaldef.h
// Martin Janson, 2000-03-09

// defines GLOBAL varaibles

//  ----- 	Crystal structure data  ---------------
// target data
int Nspecies;	// number of atomic species in lattice
struct atomStruct *atom;  // ion (atom[0]) and target atom  data
double nAtom, nAtomAmorph; //atomic denisty of target (cm^-3)
int randomtarget;

// Unit Cell
struct base3D UCbase, on2UCbase;
struct SCstruct *UCatom;
int NatUC; // Number of atoms in unit cell

// Sub Cells
struct base3D SCbase, on2SCbase;
struct int3D nSC;
int NatSC;      // total number of atoms in SC structure
int NSC; // number of subcells = nSC.x*nSC.y*nSC.z
double SCsize; // defines the number of sub cells (i.e. nSC), see "cr_subcells"
struct SCstruct *SCatom;
int 	 *SCpointer;



// -----  Primary ion data -----------------------------
struct ionStruct ion0;   // primary ion data (.E, .id, .ev0)
double dose;  // implanted dose to be simulated (cm-2)

// -----  Target surface and layer width  ---------------
struct double3D Rsurf;  // point at surface
struct double3D ensurf, e1surf, e2surf; // ensurf = surface normal
double layerWidth; //(Å)
double randsurflayer; // (Å)
double U_s; // (keV)

//  Implanted area definition ----------------------
struct double3D Rimp;
double impArea1, impArea2; // (Å, Å)

// Nuclear stopping data
int IAPfunction;

// ----- Electronic stopping power data
struct espStruct *esp;
double E_fermi;
double *locnorm;

// BCA paramters ----------------------------
double pmax, ppmax, Estop;    // (Å, Å^2, keV )
double dRmin, dRsc; // (Å)
int simcollalg, vibalg;  // flags for simulataneous collision and vibr. alg

// Rare event variables
int rareEventMult;
double rareDepth[MAXNRAREDEPTH + 1];
int NrareDepth;

// Full cascade linked list of ions
struct illStruct *ill;

// Damage model
int recoil;
double rndscatfac;

// Concentrations of primary ion,
// target interstitials and vacancies
double *Czion, *CzI, *CzV; // (cm-3)
int saveionsR3;
int saverecoilsR3;
double savedose[10];
int Nsavedose;

// Energies lost to electronic exitation Qa
// and phonons Ta (full recoil model) (energy loss in elastc collisions)
double *Qa, *Ta;  // (keV/Å/cm2)

// Mean impact parameter in each slab for each species, Fnep
int calcFnep;
double *Fnep;

// total dose (cm^-2) of backscattered and transmitted primary ions in sim.
//double totbackscatt, tottransmitt;

// Input statistics variables ---------------------
int	Nions; // number of pseudo ions in simulation
int  Ndistz; // here: z//esurf
double dz; // dz = layerWidth/nDistz  (Å)

// Output file identifyers (should later not be defined globally!!)
FILE *fionR3, *fIR3, *fVR3;



// miscellaneous
long randomz;

// -------------------------------------------------------

