// file: siimpl.c
// Martin Janson, 2000-03-10 - 2003-08-15

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "siimpl.h"
#include "mcstd.h"
#include "spstd.h"
#include "mcinit.h"
#include "mcio.h"
#include "mtrxfun.h"
#include "mjmath.h"

#include "mcstructdef.h"
#include "mcglobalextern.h"
#include "mcglobaldef.h"                                        


#define pi 3.141592654
#define EPS  1e-20
#define EPSr 1e-10
#define DRVNOTFOUND 1.0 // (1 Å added in ev0 direction if no collision
						    //  partner is found in subcell)



int main (int argc, char *argv[])
{
	char *impspecfile;
   struct ionStruct ion;
   struct illStruct *ptrsh;
	struct double3D ionR0;
	int backscatt, transmitt;
   double totbackscatt = 0.0, tottransmitt = 0.0, totsputter = 0.0;
   double dosenow = 0.0;
   double  Rp;
   double  Rpmean =  0.0 , Rpathmean = 0.0;
   double Erest, Ebacktrans = 0.0;
   double rand1, rand2;

   unsigned int totNions = 0;
   unsigned int totcol = 0, iontotcol = 0, collisions;
   unsigned int simcollev = 0;
   unsigned int totsimcol = 0;
   unsigned int nofound = 0;

	int i,ii,j;


   // read simulation file from command line ------------------
   if (argc!=2) {
     	printf("No simulation file specified!\n");
      return 1;
   }
   impspecfile=argv[1];

   // initialize simulation data ------------------------------
   initSimulation(impspecfile);

   // Open CzData files
   saveCzData(0);

   // Open file for saving 3D-positions of ions ----------------
	if (saveionsR3)
		openionR3file();

   // Open files for saving 3D-positions of self-Is and Vs
	if (saverecoilsR3)
   	openIVR3files();

	// start simulation ----------------------------------------
   for (i=0; i<Nions; i++) {

 		//  ------	Define primary pseudo ion (i) -------------
      // initial id, E, ev0, and w  from ion0 (ion.Rpath is also set to 0 here)
   	ion = ion0;

      // create set of random numbers for mj_vib algorithm
      if (vibalg == MJVIB)
      	mjvibrndnumbers(&(ion.vib));

		// find random start point for ion within impl. Area rectangle
      rand1 = randmj();
      rand2 = randmj();

      ion.r.x = Rimp.x + (impArea1*rand1*e1surf.x + impArea2*rand2*e2surf.x);
		ion.r.y = Rimp.y + (impArea1*rand1*e1surf.y + impArea2*rand2*e2surf.y);
		ion.r.z = Rimp.z + (impArea1*rand1*e1surf.z + impArea2*rand2*e2surf.z);

		// add ( epsilon * (-ensurf) ) to ion.r to ensure that ion starts
      // at the correct side of surface interface (to avoid backscatter at once)
      ion.r.x -= EPSr * ensurf.x;
		ion.r.y -= EPSr * ensurf.y;
		ion.r.z -= EPSr * ensurf.z;


      // save start 3D position of primary ion if saveionsR3
 		if (saveionsR3) ionR0 = ion.r;

      // -----  Initiate ion linked list (ill)
      ill = (struct illStruct *) malloc(sizeof(struct illStruct));

		// bottom of stack is marked with ill->ion.id = -1
      ill->ion.id = -1;

		// Add primary ion to list
		add2ill(&ion);
      //  -----


		// Start simulatation of pseudo primary ion (i) and cascade.
      do {
         // 	Read ion data from top of ion linked list
			ion = ill->ion;
         // and delete top of the stack of ill
         ptrsh = ill->p;
         free(ill);
         ill = ptrsh;

   		backscatt = 0;
   		transmitt = 0;
         collisions = 0;

			// --- Follow ion until (ion.E <= Estop)
         // ---or backscattered or transmitted
		   pseudoIonSim(&Rp, &backscatt, &transmitt, &ion, &Erest,
             &collisions, &simcollev,  &totsimcol, &nofound);

         totcol += collisions;

	      // ----  Determine Rp slab j:  Cz(Rp) = Cz[j]
	   	j = floor(Rp/dz);
         // Check that j is valid. Ion can be (just) outside the target
         // if the energy of the ion outside the target is smaller
         // than the surface energy limit so that no backscatter event
         // occurs
         if (j < 0) j = 0;
         if (j >= Ndistz) j = Ndistz - 1;


         if (ion.id == 0 )
         // Priamary ion
         {
         	++totNions; // totNions include the rare event ions
         	dosenow += ion.w;
				// Save Cz data if dosenow = savedose
            for (ii = 0; (ii < Nsavedose); ii++)

            	// Save data to files for multi dose save
            	if ( ( dosenow >= savedose[ii]) &&
                    ((dosenow - ion.w) < savedose[ii]) )
               	saveCzData(1);

         	iontotcol += collisions;

            // Update Rp and Rpath
            // NOTE only updated if the ion not backscattered.
            // this is doen to fet a correct Rp/Rpath fraction for
            // the ions left in the target. Rp and Rpath normalized
            // after mainj loop.
            if (!backscatt)  {
 		      	Rpmean += Rp * ion.w;
            	Rpathmean += ion.Rpath * ion.w;
            };

	      	if ( (!backscatt) && (!transmitt) ) {

            	Czion[j] += ion.w/(dz*1.0e-8);  // (cm-3)

					// Rest energy is trandfered to the lattice
               Ta[j] += Erest/(dz) * ion.w;  // (keV/Å/cm^2);

	      		// save start and final 3D position of ion if saveionsR3
 					if (saveionsR3)
		      		fprintf(fionR3,"%e\t%e\t%e\t\t%e\t%e\t%e\n",
                  	ionR0.x, ionR0.y, ionR0.z, ion.r.x, ion.r.y, ion.r.z);


         	}
	      	else {
            // i.e. backscattered or transmitted
            	totbackscatt += backscatt * ion.w;  //cm^-2
		     	  	tottransmitt += transmitt * ion.w;  //cm^-2
               Ebacktrans += Erest * ion.w;  // (keV/cm^2);
	      	}

         }
         else
         // i.e. when ion.id > 0 -> recoil
         // --- Update recoiling ion distriburions
         {
	      	if ( (!backscatt) && (!transmitt) ) {

	        	   CzI[j+(ion.id - 1)*Ndistz] += ion.w/(dz*1.0e-8);  // (cm-3)
					// Rest energy is trandfered to the lattice
               Ta[j] += Erest/(dz) * ion.w;  // (keV/Å/cm^2);
	      	}
	      	else {
            	totsputter += ion.w * backscatt;  //cm^-2
					// transmitted recoils are not registered
               Ebacktrans += Erest * ion.w;  // (keV/cm^2);
	      	}

 				if (saverecoilsR3)
		      	fprintf(fIR3,"%d\t%e\t%e\t%e\n",
               					ion.id, ion.r.x, ion.r.y, ion.r.z);
         } // ---------------------------------

		} while (ill->ion.id !=-1);

      free(ill);
   }; // end of: for (i=0; i < Nions ...


   // Simulation completed -----------------------------------------

	// Normalize Rp and R

   Rpmean /= (dose - totbackscatt);
   Rpathmean /= (dose - totbackscatt);


   // Close R3 files
   if (saveionsR3) fclose(fionR3);
   if (saverecoilsR3){
   	 fclose(fIR3);
   	 fclose(fVR3);
   }

   // Write CzData files
   saveCzData(1);

   // Close CzData files
   saveCzData(2);

   // Write simulation data to log file
   writelogfile(Rpmean, Rpathmean,
   				totbackscatt, tottransmitt,	totsputter, Ebacktrans,
              	totcol, iontotcol,  simcollev, totsimcol, nofound, totNions);



 	// Exit simulation program
   freeAlocatedMemory();
	printf("\nSimulation completed successfully!\nSiiMPL - (c) M.S.Janson, 2003\n\n");
   return 0;
} // end of: main
// -----------------------------------------------------------------------------



void pseudoIonSim(double *Rp,
						int *backscatt, int *transmitt,
					   struct ionStruct *ion, double *Erest,
                  unsigned int *collisions,
                  unsigned int *simcollev, unsigned int *totsimcol,
                  unsigned int *nofound)

{
	struct double3D R2UC, RcUC, RcUClast, ev1, ev2;
	struct double3D imp, tmpdbl3D;
   struct double3D tgatomlist[MAXSCATOMS];
   int	 Ntgatom, lastSCi, amorphous, damagescatter;
   struct int3D UCi, lastUCi;
   double Se, Q, p, p0, b, T, Tlatt, THETA, M1, M2, Prnd;
	double dRv0, dRv0last, tmpdbl;
   double SeLo, SeHi, f, p_esp;
   struct ionStruct recoilion;
	struct scStruct sclist[MAXNSIMCOLL + 1]; // "+1" since first bc is stored
   												     // at the end of sclist (see below)

	int rareDepthPoint = 0;

	#define NUD_DR 50
   double dR0, dRloc, dRloc_ud;
	int i_ud_dR = 0;

   int id2, Z1, Z2, tgatomi, SCp1, Nsimcoll;
	int i, zi;


   lastSCi = -1;

   M1 = atom[ion->id].M;
   Z1 = atom[ion->id].Z;

   dRmin = 0.0;

   // Damage variables
   dR0 =  (1.0/(nAtom*1e-24)/(pi*sqr(pmax)));
   dRloc = dR0;
   dRloc_ud = 0.0;

   // 	Calculate initial projected range of ion
   //  i.e. Rp != 0 for recoils and rare event primary ions
   // Rp = (-ensurf) dot (Rion - Rsurf)
   *Rp =   (-ensurf.x) * (ion->r.x - Rsurf.x)
         + (-ensurf.y) * (ion->r.y - Rsurf.y)
         + (-ensurf.z) * (ion->r.z - Rsurf.z);

   zi = floor(*Rp/dz);
   // Ckeck of zi (for recoils)
   if (zi<0) zi = 0;
   if (zi>=Ndistz) zi = (Ndistz - 1);

   // Make sure that the ion is not backscattered ...
   if ( (*Rp < 0.0) &&
          (( (ensurf.x) * (ion->ev0.x) +
             (ensurf.y) * (ion->ev0.y) +
             (ensurf.z) * (ion->ev0.z)   ) > 0.0) )
   {
   	*backscatt = 1;
   }
   // ... or transmitted from the beginning. This could happen for recoils only
   else if (*Rp >= layerWidth) *transmitt = 1;

   // Determine rareDepthPoint, i.e. between which rareDepths that
   // the ion is possitioned, != 0 for rareEvent ions
   // For primary ion only
   if ( (ion->id) == 0) {
   	rareDepthPoint = 0;
      while (*Rp > rareDepth[rareDepthPoint]) rareDepthPoint++;
   }

   // ----------- Start main ion loop ------------------------------------
   while ( (ion->E >= Estop) && (! *backscatt) && (! *transmitt) )
   {

   	Nsimcoll = 0;

   	// Calculate ion.UCr0, ion.SCi, ion.rUC and UCi for ion.r
   	r2ucsc(ion, &UCi);

		// --- Find collision partner in random or cryst. target

      // Probability for scatter at damage,
      //    Prnd = rndscatfac * sum(CzI(i)) / nAtom * <dRloc>/dR0;
      Prnd = 0;
      for (i=0; i<Nspecies; i++)
      	Prnd += CzI[zi + i * Ndistz];
      Prnd *= (rndscatfac / nAtom);

      // Correction term to P
      Prnd *=  dRloc/dR0;

      amorphous = 0;
      damagescatter = 0;
      // Determine if target should be random
      if ( (randomtarget) || (*Rp < randsurflayer) ) {
          // Amorphous target  => n = nAtomAmorph
   		findrandomtg(&R2UC, &RcUC, &dRv0, &p, &tgatomi,
         	/*in*/  ion, &nAtomAmorph);
         amorphous = 1;
      }
      else if (randmj() < Prnd ) {
      	// Scatter at damage in crystal => n = nAtom
   		findrandomtg(&R2UC, &RcUC, &dRv0, &p, &tgatomi,
         	/*in*/  ion, &nAtom);
         damagescatter = 1;
      }
      else  {
      // Crystalline target
      	if (!simcollalg)
         {
         	switch (vibalg)
         	{
         		case ROBINSONVIB:
		         	findbc_robvib(&R2UC, &RcUC, &dRv0, &p, &tgatomi, /*in*/ ion);
         			break;
         		case MJVIB:

               	// Update tgatomlist if SC (SCi and UCi) have changed
               	if ( !( (ion->SCi == lastSCi) &&
                          (UCi.x == lastUCi.x) &&
                          (UCi.y == lastUCi.y) &&
                          (UCi.z == lastUCi.z)  )  ) {

                  	updatetgatoml_mjvib(tgatomlist, &Ntgatom, &SCp1, /*in*/ ion);
                     lastSCi = ion->SCi;
                     lastUCi = UCi;
                  };

 		         	findbc_mjvib(&R2UC, &RcUC, &dRv0, &p, &tgatomi,
                  /*in*/ tgatomlist, &Ntgatom, ion);

      		     // NOTE when tgatomlist is used, tgatomi is returned as the
                  // number in that list. To relate to the index in the SCatom
                  // list, 'SCp1', (from updatetgatom_mjvib) is added to tgatomi
                   if (tgatomi != -1)
                  	tgatomi += SCp1;


					break;

         		default:
         			printf("%d is not a valid vibalg!\n",vibalg);
         			exit(1);
            }
         }
         else // i.e. simultaneous collision algorithm
         {
         	switch (vibalg)
         	{
         		case ROBINSONVIB:
	         		findbc_sc_robvib(&R2UC,  &RcUC,  &dRv0,  &p,  &tgatomi,
			                        sclist, &Nsimcoll, /*in*/ ion);
               break;
         		case MJVIB:
               	// Note that findbc_sc_mjvib do not (yet) use the tgatom list
                  // and therefore give the indexing relative SCatom[i]
	         		findbc_sc_mjvib(&R2UC,  &RcUC,  &dRv0,  &p,  &tgatomi,
			                        sclist, &Nsimcoll, /*in*/ ion);
               break;
         		default:
         			printf("%d is not a valid vibalg with sim. coll. alg.!\n"
                  ,vibalg);
         			exit(1);
            }
				if (Nsimcoll != 0)
            	++(*simcollev);
            	(*totsimcol) += (Nsimcoll + 1);
         }
      } // end if: Crystalline target

      // Add first collision to "sim.coll.list" sclist
      // i.e. if only binary collision, Nsimcoll = 1
      // and collision data are stored in sclist[0]
      sclist[Nsimcoll].tgatomi = tgatomi;
      sclist[Nsimcoll].RcUC = RcUC;
      sclist[Nsimcoll].R2UC = R2UC;
      sclist[Nsimcoll].p = p;
      sclist[Nsimcoll].dRv0 = dRv0;
      ++Nsimcoll;


      dRv0last = 0.0;

      Q = 0.0;
		Tlatt = 0.0;
      imp.x = 0.0;
      imp.y = 0.0;
      imp.z = 0.0;


       // ---  Total avarage Electronic stopping Se  -----

       //Se_lo = A1 * E^(p + f*p2),   from SiC stopping paper
       if (ion->E > 1.0)
       	f = exp(-3.0* pow((log(E_fermi * atom[ion->id].M )/log(ion->E)), 4.0));
       else
       	f = 0.0;

       p_esp = esp[ion->id].p + f * esp[ion->id].p2;
       SeLo = esp[ion->id].A1 * pow(ion->E, p_esp);

       // Se_hi = A3 / E * log(1 + A4/E + A5*E)
       // (NOTE (E/1000) included in A4-A5 parameters) 
       SeHi = esp[ion->id].A3/ion->E *
	       log( 1.0 + esp[ion->id].A4/ion->E + esp[ion->id].A5*ion->E );

       // Se = (Se_lo * Se_hi)/(Se_lo + Se_hi)
       if ((SeHi > EPS) & (SeLo > EPS) )
       	Se = 1.0/( 1.0/SeLo + 1.0/SeHi);
       else
       	Se = SeLo;
       // -----------------------------------------------


      // ----- Local energy losses, Q, T, (imp) ,(loop for sim. coll)

   	// Make ceratain that collision partner is found in Sub-Cell
    	if ( tgatomi != -1 )
      {
      	// Note that collisions not (necessary) are taken in right order.
      	for (i=0; i < Nsimcoll; i++)
      	{


      		++(*collisions);

	      	// Update target data for next (simultaneous) collision
	      	tgatomi = sclist[i].tgatomi;
	      	R2UC 	  = sclist[i].R2UC;
	      	RcUC 	  = sclist[i].RcUC;
	      	dRv0 	  = sclist[i].dRv0;
	      	p 	 	  = sclist[i].p;

				// Find last collision values of dRv0 and RcUC
            // i.e. the collison furtherst away (of all simultaneous collisions)
            // relative to the current ion position
				if (dRv0>dRv0last)
         	{
         		dRv0last = dRv0;
					RcUClast = RcUC;
         	}

            // target atom data
      		id2 = SCatom[tgatomi].id;
				Z2 = atom[id2].Z;
				M2 = atom[id2].M;

         	b = p/aZBL(Z1, Z2);

            if (calcFnep) {

            	if (!(amorphous || damagescatter)) {
            		// Determine impactparameter to equilibrium (non-vibrating)
               	// target atom

			   		// tmpdbl3D (rIonTg) = vector going from ion to target atom i
               	tmpdbl3D = sub3D( &SCatom[tgatomi].r, &ion->rUC);

      				// tmpdbl (dR) = projection of rIonTg on velocity axis
               	//	= (rIonTg)dot(ev0)
      				tmpdbl = dot3D( &tmpdbl3D, &ion->ev0);

	      			// calculate impact parameter, p0
         			p0 = 	sqrt(
               		sqr(tmpdbl3D.x - tmpdbl*(ion->ev0.x) ) +
                  	sqr(tmpdbl3D.y - tmpdbl*(ion->ev0.y) ) +
                  	sqr(tmpdbl3D.z - tmpdbl*(ion->ev0.z) ) );

               }
               else
               // collision with 'random' atom
               	p0 = p;

               // Random mean area per atom
               if (amorphous)
                  tmpdbl = pow( 1.0e24/nAtomAmorph, 0.66667);  // Å^2
               else
               	tmpdbl = pow( 1.0e24/nAtom, 0.66667);  // Å^2

            	// Update nuclear encounter probability distribution
            	// (Fnep is normalized aginst the ion pseudo dose)
            	// F = 1/(2*pi*u1^2)*exp(-p^2/(2*u1^2)) * A
               //	  =  3/(2*pi*u3^2)*exp(- 3/2*(p/u3)^2) * A
               if (atom[id2].u3 > 0.0)
            		Fnep[zi+(id2 - 1)*Ndistz] +=
            			3.0/( 2.0*pi*sqr(atom[id2].u3) )
               		* exp(- 3.0/2.0 * sqr(p0/atom[id2].u3) )
                     * tmpdbl* ion->w/dose/dz;
            };

				//  ---  Local Electronic stopping for crystal target
            if (!amorphous)
	     			Q += esp[ion->id].lf * Se
            		* locnorm[ion->id * Nspecies + (id2 - 1)]
               	* exp(-b*esp[ion->id].s) ;


         	// ----	Nuclear stopping and scattering ----

         	// Determine scattering angle THETA in CM-units
         	switch (IAPfunction)
         	{
         		case MAGICFORMULA:
         			//	magicFormula gives:
         			// ans = sin(THETA/2)^2 = [1-cos(THETA)]/2;
                  THETA = acos( 1.0 - ( 2.0
         				* magicFormula(keV2epsl(ion->E, Z1, M1, Z2, M2), b)) );
         			break;
         		case GAUSSMEHLERZBL:
         			THETA = ThetaGMzbl( keV2epsl(ion->E, Z1, M1, Z2, M2), b) ;
         			break;
         		default:
         			printf("%d is not a valid IAPfunction!\n",IAPfunction);
         			exit(1);
         	}

         	// 	Calculate scatter directions (ion ev1, target ev2)
         	// and elastic energy loss T in binary collision
         	bcScatter(&ev1, &T, &ev2,
         			/*in*/ ion, &R2UC, M2, &RcUC, THETA);

         	// Add recoil impulse to total impulse
         	tmpdbl = sqrt( 2.0 * M2 * T);
         	imp.x += ev2.x * tmpdbl;
         	imp.y += ev2.y * tmpdbl;
         	imp.z += ev2.z * tmpdbl;

	      	// ------------- Full cascade damage model
      		if (recoil)
         	// NOTE that the total recoiling energy is larger than
         	// the energy loss by the ion in a simulataneous collision
         	// The recoiling energies should be scaled in a later version of siimpl
      		{
         		// Stack recoiling atoms and create vacancy if T >= Edispl
         		if (T >= atom[id2].Ed)
         		{
         			// Position of recoiling atom = rSC +  (EPSr * ev2)
            		// (to avoid self-collision)
               	// NOTE that the start position of the recoiling atoms must be
               	// treated different in the MJVIB mode where collision partners
               	// are picked WITH vibration added.
               	switch (vibalg)
         			{
         			case MJVIB:
		           		recoilion.r.x = R2UC.x + (ion->UCr0.x)
      	      			+ ( EPSr * ev2.x );
         	   		recoilion.r.y = R2UC.y + (ion->UCr0.y)
            				+ ( EPSr * ev2.y );
               		recoilion.r.z = R2UC.z + (ion->UCr0.z)
            				+ ( EPSr * ev2.z );
        				break;
         			default: // i.e. for ROBINSONVIB and PROJVIB
 		           		recoilion.r.x = SCatom[tgatomi].r.x + (ion->UCr0.x)
      	      			+ ( EPSr * ev2.x );
         	   		recoilion.r.y = SCatom[tgatomi].r.y + (ion->UCr0.y)
            				+ ( EPSr * ev2.y );
               		recoilion.r.z = SCatom[tgatomi].r.z + (ion->UCr0.z)
            				+ ( EPSr * ev2.z );
            		}

            		recoilion.id = id2;
            		recoilion.E = T;
            		recoilion.ev0 = ev2;
            		recoilion.w = ion->w;

            		// Add recoiling ion to ion linked list ill
            		add2ill(&recoilion);

            		// Update vacancy distribution
            		CzV[zi+(id2 - 1)*Ndistz] += ion->w/(dz*1.0e-8);  // (cm-3)

            		// Save R3 postion of vacancies
            		if (saverecoilsR3)
            			fprintf(fVR3,"%d\t%e\t%e\t%e\n",
                     recoilion.id, recoilion.r.x, recoilion.r.y, recoilion.r.z);

          		}
            	else
               // i.e. if energy is not large enough to create recoil
               Tlatt += T;

      		}// end of: Full cascade damage model

         	else
         	// Kinchin-Pease damage model
   	   	{
      			// Update interstitial and vacancy distributions
               tmpdbl = 0.4 * T /(atom[id2].Ed) * (ion->w)/(dz*1.0e-8); //(cm-3)
      			CzI[zi+(id2 - 1)*Ndistz] += tmpdbl;
      			CzV[zi+(id2 - 1)*Ndistz] += tmpdbl;

 	      		Tlatt += T;
				}

         } 	// end of: local energy loss loop, Nsimcoll
      } // end of if (!tgatomi)
      else
      // i.e. if no collision partner is found in sub-cell
      {
   			++(*nofound);
      		dRv0last = DRVNOTFOUND;
       		T = 0.0;
            RcUClast.x = ion->rUC.x + dRv0last*ion->ev0.x;
            RcUClast.y = ion->rUC.y + dRv0last*ion->ev0.y;
            RcUClast.z = ion->rUC.z + dRv0last*ion->ev0.z;
      }

   	// Add non-local electronic stopping power to total Q or
      // ,if amorpous target, total Q
      if (!amorphous)
      	// dQ = (1 - lf)*Se*dR
  			Q += (1.0 - esp[ion->id].lf) * Se
     			* nAtom*1e-24 * dRv0last;
      else
  			Q = Se * nAtomAmorph*1e-24 * dRv0last;

		// --- Calculate NSP energy loss T
      //     and ion scattering direction ev0 from imp

      // Subtract total target impulse from ion impulse, imp0
      tmpdbl = sqrt( 2.0 * M1 * ion->E);  // tmpdbl = abs(imp0)
      imp.x = (tmpdbl * ion->ev0.x) - imp.x;
      imp.y = (tmpdbl * ion->ev0.y) - imp.y;
      imp.z = (tmpdbl * ion->ev0.z) - imp.z;

      // Calculate new ion direction ev0
      tmpdbl = norm3D(&imp);  // tmpdbl = abs(imp0 - sum(imp) )
      ion->ev0 = divk3D(&imp, tmpdbl);

      // Calculate ion->E = |p1|^2/(2*M1)
      // Equal to "ion->E -= T " in a true binary collision
     	ion->E = sqr(tmpdbl)/(2.0 * M1 );

		if ( ion->E < 0.0 )
      // This can happen in the simultaneous collision alg.
   	{
      	Tlatt += (ion->E);
         Q = 0.0;
      	ion->E = 0.0;
      }
      // ----------------------------------

      // --- Subtract total Q from ion->E
      ion->E -= Q;

		if ( ion->E < 0.0 )
      // This can happen for low ion->E and small p (with local ESP)
   	{
      	Q += (ion->E);
      	ion->E = 0.0;
      }

      // --- Update ion position
      ion->r = add3D(&RcUClast, &ion->UCr0 );  // r = rUC + UCr0

    	ion->Rpath += dRv0last;

      // 	Calculate projected range of ion
      // Rp = (-ensurf) dot (Rion - Rsurf)
      *Rp =   (-ensurf.x) * (ion->r.x - Rsurf.x) +
      		  (-ensurf.y) * (ion->r.y - Rsurf.y) +
              (-ensurf.z) * (ion->r.z - Rsurf.z);
      zi = floor(*Rp/dz);
      // Check of zi
      if (zi<0) zi = 0;
      if (zi>=Ndistz) zi = (Ndistz - 1);

		// Add to local  dRloc_ud (for calculation of local dR)
      dRloc_ud += dRv0last/NUD_DR;
      i_ud_dR++;

      // Update dRloc every NUD_DR collision
      if (i_ud_dR == NUD_DR ) {
      	dRloc = dRloc_ud;
         dRloc_ud = 0.0;
         i_ud_dR = 0;
      };


      // update energy loss distributions
      Qa[zi] += Q / dz * (ion->w);  // (keV/Å/cm^2)
      Ta[zi] += Tlatt / dz * (ion->w);  // (keV/Å/cm^2)


      // Check if ion is backscattered, (or sputtered) ...
      //    the ion must be directed out from the surface for
      //		a true backscatter event. Due to (ROBINSON) vibration algorithm.
      if ( (*Rp < 0.0) &&
           (( (ensurf.x) * (ion->ev0.x) +
              (ensurf.y) * (ion->ev0.y) +
              (ensurf.z) * (ion->ev0.z)   ) > 0.0) )
  			*backscatt = 1;
		// ... or transmitted
   	else if (*Rp >= layerWidth)
      	*transmitt = 1;


      // Determine largest dRmin of all collisions in  (1 for bc)
      // (dRmin is used in the Robinson vibration algorithm as well
      //  as in all sim. coll. algs. to avoid self collision)
      dRmin = DRMINMIN;
      for (i = 0; i<Nsimcoll; ++i)
      {
         if (vibalg == ROBINSONVIB) {
         // In the Robinsin alg. search for next collision is done
         // by lattice position (without vibration)
         tmpdbl =
         	(SCatom[sclist[i].tgatomi].r.x - sclist[i].RcUC.x ) * (ion->ev0.x) +
            (SCatom[sclist[i].tgatomi].r.y - sclist[i].RcUC.y ) * (ion->ev0.y) +
            (SCatom[sclist[i].tgatomi].r.z - sclist[i].RcUC.z ) * (ion->ev0.z);
	         if (tmpdbl > dRmin) dRmin = tmpdbl;
         }
         else if ( (vibalg == MJVIB) && (simcollalg)){
         // In the mj vib sim coll. alg. search for next collision is done
         // by atom position including (reproducable) vibration.
         tmpdbl =
         	(sclist[i].R2UC.x - sclist[i].RcUC.x ) * (ion->ev0.x) +
            (sclist[i].R2UC.y - sclist[i].RcUC.y ) * (ion->ev0.y) +
            (sclist[i].R2UC.z - sclist[i].RcUC.z ) * (ion->ev0.z);
	         if (tmpdbl > dRmin) dRmin = tmpdbl;
         }
      }

      // Check for rare depth event
      // Rare depth only evaluated for primary ion
      if ( (ion->id) == 0) {

      	if  (*Rp > rareDepth[rareDepthPoint]) {
         	++rareDepthPoint;
   	   	rareEvent(ion);
      	}
      }

	} //  end: while ( (ion->E >= Estop) && (! *backscatt) && (! *transmitt) )

   
	// Check that backscattered ion has enough energy to pass the surface barrier
   if (*backscatt && (ion->E < U_s) )
   	*backscatt = 0;

	*Erest = ion->E; // (keV)

} // end of: pseudoIonSim  -----------------------------------------



void add2ill(struct ionStruct *ion)
{
   struct illStruct *ptrsh;

   ptrsh = ill;
   if ((ill = (struct illStruct *) malloc(sizeof(struct illStruct)) )== NULL)
   {
   	printf("Out of memory in add2ill function!\n");
      exit(1);
   }
   ill->p = ptrsh;
   ill->ion = *ion;
}
//  -----------------------------------------------------------





