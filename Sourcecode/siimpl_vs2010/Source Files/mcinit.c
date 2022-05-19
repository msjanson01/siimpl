// file: mcinit.c
// Martin Janson, 2000-03-20

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mcinit.h"
#include "mcio.h"
#include "mtrxfun.h"
#include "spstd.h"
#include "mjmath.h"
#include "subcells.h"

#include "mcstructdef.h"
#include "mcglobalextern.h"

#define pi 3.141592654


int initSimulation(char *filename)
// Set default simulation values ans
// initialize simulation defined in "filename"
{
   struct double3D a3D;
   int i,j,count;
   double a, Z2eff, I, C;


   // default Ndistz and Nspecies value
   Ndistz = 100;

   Nspecies = 0;
   NatSC = 0;
	NSC = 0;

   // -----   Look for structure size data in "filename"
	structureSize(filename, &NatUC, &Nspecies, &Ndistz);

   if ( (NatUC == 0)||(Nspecies==0) ){
   	printf("Error in 'UCatom' or 'Z2' definitions!\n");
      printf("NatUC = %d,  Nspecies = %d\n", NatUC, Nspecies);
      exit(1);
	}


	// ----- 	Allocate memory for Global variables ----------------
   // structure variables
   UCatom = (struct SCstruct *) malloc(NatUC*sizeof(struct SCstruct));

   atom = (struct atomStruct *) malloc((Nspecies+1)*sizeof(struct atomStruct));

	//  electronic stopping power
   esp = (struct espStruct *) malloc((Nspecies+1)*sizeof(struct espStruct));
   locnorm = (double *) malloc((Nspecies+1)*(Nspecies)*sizeof(double));

	//	ion, damage, energy, p_mean distribution variables
   Czion = (double *) malloc(Ndistz*sizeof(double));
   CzI = (double *) malloc(Ndistz*Nspecies*sizeof(double));
   CzV = (double *) malloc(Ndistz*Nspecies*sizeof(double));

   Qa = (double *) malloc(Ndistz*sizeof(double));
   Ta = (double *) malloc(Ndistz*sizeof(double));

   Fnep = (double *) malloc(Ndistz*Nspecies*sizeof(double));


   // -----    Read simulation data from file "filename" ----------
	// This first reading of 'filename' is only done so that information
   // of the crystal is read into 'atom' to be used in the calculation
   // of the ESP  default values below and in next readSimData (if esp
   // data is defined before the structure data in 'filename')
   readSimData(filename);

   // Calculate relative abundance x of target atom i
	for (i = 1; i <= Nspecies; i++) {
   	count = 0;
   	for (j = 0; j < NatUC; j++)
      	if  (UCatom[j].id == i ) ++count;

      atom[i].x = (double) count/NatUC;
   }

   // Determine atomic density in (crystalline) target
   nAtom = NatUC / volumebase3D(UCbase) * 1e24; // (atoms/cm^3)

	// -----  Set defualt values -------------------------------------

   // Target atom data
   for (i= 1; i<= Nspecies; i++){
   	atom[i].Ed = 15; // eV
   	atom[i].u3 = 0.0; // Å
   }

   for (i=0; i< NatUC; i++){
   	UCatom[i].id = 1;
   }

   // Atomic density for amorpous material
	nAtomAmorph = nAtom;

   // Sub Cell size, used for determineing the nuber of Sub Cells in cr_subcells
   SCsize = 2.0;

   // primary ion (ion0) data
   ion0.id = 0; // primary ion always has '.id=0' !
   ion0.ev0.x=0.0; ion0.ev0.y=0.0; ion0.ev0.z=1.0;
   ion0.Rpath = 0.0;

	// target surface and layer width
	Rsurf.x=0.0;  Rsurf.y= 0.0;  Rsurf.z= 0.0;
	ensurf.x=0.0; ensurf.y= 0.0; ensurf.z=-1.0;
	layerWidth = 10000.0; // (Å)
   U_s = 5e-3; // (kev)

   randsurflayer = 0.0; //(Å)

	randomtarget = 0;

	//  Implanted area definition
	Rimp = Rsurf;
	impArea1 = 1000.0; // (Å)
   impArea2 = 1000.0; // (Å)

   // Nuclear stopping data
	IAPfunction = MAGICFORMULA;

   // Rare event models
   rareEventMult = 10;
   // NOTE that rareDepth[NrareDepth] is set to Inf at end of init
	NrareDepth = 0;


   // ----- Electronic stopping  power
   E_fermi = 25.0; // (keV/amu)


  	// Determine effective Z2 for BB stopping
  	Z2eff = 0.0;
  	for (j=1; j<=Nspecies; j++)
   	Z2eff += atom[j].x *atom[j].Z;

   // Determine I from Z2 (Z2eff) ,according to Biersack
   if (Z2eff < 13.0)
   	I = Z2eff * (12.0 + 7.0/Z2eff) * 1.6022e-012; // (erg)
	else
   	I = Z2eff * (9.76 + 68.5 * pow(Z2eff, -1.19)) * 1.6022e-012; // (erg);


   for (i=0; i<=Nspecies; i++) {

   	// --- Low energy parameters
   	esp[i].p = 	0.5;
   	esp[i].p2 = 0.0;
	   esp[i].s = 	0.3;
      esp[i].lf = 1.0;

      // Low velocity Se according to LS and Braggs rule
      esp[i].A1 = 0.0;
      for (j=1; j<=Nspecies; j++) {
      	esp[i].A1 +=
         		(atom[j].x) * SeLS(1.0, atom[i].Z, atom[i].M, atom[j].Z);
      }

      // -- High energy stopping papermeters according to Biersack's formula
      // (NOTE (E/1000) included in A4-A5 parameters)
 		// A3 = a3 * Z1^2*Z2*M1 (keV*Å^2)
 	  	esp[i].A3 = 2.38 * sqr(atom[i].Z) * Z2eff * atom[i].M;

      // Determine C = C(Z1, Z2), Biersack
	   if (atom[i].Z < 3)
   		C = 100.0 * atom[i].Z / Z2eff;
		else
   		C = 5.0;

    	esp[i].A4 = 2.84e11 * C * I * atom[i].M;
      esp[i].A5 = 3.52e-12 / (I * atom[i].M);
   }

	// BCA parameters
   pmax = 2.0;   //(Å)
   Estop = 2e-3; //(keV)
   simcollalg = 0;
   vibalg = MJVIB;
   dRsc = 0.2; // Å

   // Damage model
   recoil = 0;
   rndscatfac = 0.0;

	// Input statistics variables
	Nions = 10000;

   // Set distibution and sim. result variables to 0 (unnecessary?)
   for (i=0; i<Ndistz; i++) {
   	Czion[i]=0.0;
      Qa[i]=0.0;
      Ta[i]=0.0;
   }
   for (i=0; i< (Ndistz * Nspecies); i++){
        	CzI[i]=0.0;
        	CzV[i]=0.0;
         Fnep[i] = 0.0;
   }


   calcFnep = 0;
	saveionsR3 = 0;
   saverecoilsR3 = 0;

   Nsavedose = 0;

   // Miscellaneous
   randomz = 314159L;
	// ------------------------------------------------------------


   // -----    Read simulation data from file "filename" ----------
   readSimData(filename);


   // Recalculate Unit Cell coordinates from UCbase to ortho normal base
   for (i=0; i< NatUC; i++)
   	UCatom[i].r = vecttransf( UCatom[i].r, UCbase);

   // on2UCbase
	on2UCbase = invbase(UCbase);

   // Calculate subcell structure
   cr_subcells();

   // Calculate squred pmax
   ppmax = pmax*pmax;

   // --- Electronic Stopping Power  ------

   // locnorm, Normalizing constant for local esp model
   // one value is calculated for each ion (primary + recoils) and atom comb.
   for (i=0; i<=Nspecies; i++) // i = id of ion

   	for (j=1; j<=Nspecies; j++) // i = id of target atom
      {
			a = aZBL(atom[i].Z, atom[j].Z); // (Å)

        	locnorm[i * Nspecies + (j-1)]  =
         	sqr(esp[i].s) / ( 2.0*pi*sqr( a ) )
        		/ ( 1.0 - ( (1.0 + esp[i].s*pmax/a ) * exp(-esp[i].s * pmax/a) ) );
   	}
   // ---------------

   // dose carried by each pseudo ion (.w)
   ion0.w = dose/Nions;

   // find vectors e1surf, e2surf in surface plane (perp. to ensurf)
   //  a3D = vector that very unlikely will be parrallel to ensurf
   a3D.x = 1.0; a3D.y = pi; a3D.z = 1905.0/1632;

   e1surf = cross3D(&ensurf, &a3D);
   e1surf = divk3D(&e1surf, norm3D(&e1surf));
   e2surf = cross3D(&e1surf, &ensurf);

   // statistics
   dz = layerWidth/Ndistz;

   // rare event
   rareDepth[NrareDepth] = INF;


   // Warning situations
   if (ion0.E < Estop) {
   	printf("ion0.E < Estop!\n");
      exit(1);
   }

   // SaverecoilsR3 can only be set if recoils are simulated (recoil == 1)
   if ( (!recoil) && saverecoilsR3 ) {
   	printf("WARNING, 'saveIVR3' is activated");
      printf(" without damage cascade simulation mode!\n");
   	printf("Simulation terminated\n\n");
 		exit(1);
   };

   if (recoil && ( Estop > U_s)){
   	printf("WARNING, 'Estop' is larger than 'surfaceEnergy'!\n");
      printf("This will lead to an erroneous sputtering fraction\n");
      printf("'Estop' is set to 0.9*'surfaceEnergy'\n\n");
      Estop = 0.9 * U_s;
   };


   return 1;
} // initSimulation -----------------------------------------------------------



void freeAlocatedMemory(void)
{
	// Allocated in subcells.c
	free(SCatom);
   free(SCpointer);

   // Allocated in mcinit.c
   free(UCatom);
	free(atom);

   free(esp);
   free(locnorm);

   free(Czion);
   free(CzI);
   free(CzV);

   free(Qa);
   free(Ta);

   free(Fnep);

}


