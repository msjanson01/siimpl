// file: mcstd.c
// Martin Janson, 2000-03-01

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mcstd.h"
#include "mtrxfun.h"
#include "mjmath.h"
#include "siimpl.h"

#include "mcstructdef.h"
#include "mcglobalextern.h"

#define pi 3.141592654
#define EPS 1e-15

// standard binary collision functions -----------------


// Calculate UC and SC indices from R(x,y,z)
void r2ucsc(struct ionStruct *ion, struct int3D *UCi)
{
	struct double3D UC,SC;

  // Calculate rUC, r0UC, SCi
  SC = floorvecttransf(ion->r,on2SCbase); // 9 mult + 6 add

	UC.x = floor(SC.x/nSC.x);          //  3 mult
	UC.y = floor(SC.y/nSC.y);
	UC.z = floor(SC.z/nSC.z);

  // Somewhat stupid strucure here, since vecttrans needs double3D
  // Change later if 'tgatoml' is kept 
  UCi->x = (int) UC.x;
  UCi->y = (int) UC.y;
	UCi->z = (int) UC.z;

   SC.x -= UC.x*nSC.x;             //  3 mult + 3 add
   SC.y -= UC.y*nSC.y;
   SC.z -= UC.z*nSC.z;

   ion->UCr0 = vecttransf(UC,UCbase);     // 9 mult + 6 add
   ion->rUC = sub3D(&ion->r, &ion->UCr0);    // 3 add

   ion->SCi = (int) ( SC.x + SC.y*nSC.x + SC.z*nSC.x*nSC.y );     //  3 mult + 2 add

	// TEST, valid SCi: 0 -> NSC-1
   // removed from code after test period
	if ( ((ion->SCi) < 0) || ((ion->SCi) >= NSC) )
   {
		printf("Fatal ERROR in function: r2ucsc(&ion)\n");
   	printf("The program has calculated an invalid SC index\n");
   	printf("SCi = %d, NSC = %d\n",ion->SCi,NSC);
    exit(1);
   } // ---

} // end of: r2ucsc -----------



void updatetgatoml_mjvib(struct double3D *tgatomlist, int *Ntgatom, int *SCp1,
             const struct ionStruct *ion)
{
	int i, SCp2;
   double u3, phi, Theta;

	double R;

   union rndR_union {
   	double R;           // 64 bit double
      unsigned int i[2];  // i[0] == lo (32) byte of R, i[1] == hi byte of R
	};

   union rndR_union rndR[4];


   // Loop all atoms in SC to find collision partner
	*SCp1 = SCpointer[ion->SCi];
   SCp2 = SCpointer[(ion->SCi) + 1];

   *Ntgatom = SCp2 - *SCp1;

   for (i = *SCp1; (i < SCp2) ; ++i)   {

     	// --- Add vibration to target atom i at position

   	// Determine four quasi random numbers (here between 1.0 and (2 - 2^-X))
      // from R and random ion.vib numbers.
      // R = |SCatom.r + ion->UCr0|, RR = |R|^2
      R =  sqrt( sqr(SCatom[i].r.x + ion->UCr0.x) +
           		  sqr(SCatom[i].r.y + ion->UCr0.y) +
                 sqr(SCatom[i].r.z + ion->UCr0.z) + EPS );

	   // Construct a positive double rndR with exponent 0 and with
      // the X lowest mantissa bits of R shifted to the X highest
      // bits of rndR. Rlo should then be shiftged (32-1-11 - X) bits
      // and stored in rndR_hi. 0x3ff in rndR_hi unsures positive, 0 exponent
      // double. 10 lowest bits give a shift of 10.
      #define BITSHIFT 10

   	rndR[0].R = (ion->vib.A) * R;
	   rndR[0].i[1] =  ( 0x3ff00000 | ((rndR[0].i[0]<<BITSHIFT) & 0x000fffff) );

   	rndR[1].R = (ion->vib.k.x) * R;
	   rndR[1].i[1] =  ( 0x3ff00000 | ((rndR[1].i[0]<<BITSHIFT) & 0x000fffff) );

   	rndR[2].R = (ion->vib.k.y) * R;
	   rndR[2].i[1] =  ( 0x3ff00000 | ((rndR[2].i[0]<<BITSHIFT) & 0x000fffff) );

   	rndR[3].R = (ion->vib.k.z) * R;
	   rndR[3].i[1] =  ( 0x3ff00000 | ((rndR[3].i[0]<<BITSHIFT) & 0x000fffff) );


      // Calculate normal distributed vibration amplitede u3 (see Xnorm3D)
      // u3 = u30 * sqrt(-2*log(rnd)) * cos(2*pi*rnd)
   	u3 = atom[SCatom[i].id].u3 *
      	  	sqrt( -2.0 * log( rndR[0].R - 1.0 + EPS ) ) *
   		  	cos( 2.0 * pi * (rndR[1].R - 1.0) );

		// 3D random normal vector, evrand, defined by angles: phi and Theta
      // phi = 2*pi * rnd
		phi = 2.0 * pi * (rndR[2].R - 1.0);
      // Theta = acos( 2*rnd - 1)
		Theta = acos( 2.0 * rndR[3].R  - 3.0 );

      // Position of atom including vibration (relative UC)
   	tgatomlist[i - *SCp1].x = SCatom[i].r.x + u3 * sin(Theta)*cos(phi);
   	tgatomlist[i - *SCp1].y = SCatom[i].r.y + u3 * sin(Theta)*sin(phi);
   	tgatomlist[i - *SCp1].z = SCatom[i].r.z + u3 * cos(Theta);

	} // end of SC list
};// end of: updatetgatoml_mjvib




//	 ---	findbc_mjvib
// Vibration = vib(R) added to each target atom before search.
// in data:  *tgatomlist, *Ntgstom, *ion
// out data: *R2UC, position of target atom relative UC
//           *RcUC, point of collision relative UC
//				 *p, imapct parameter
//				 *dRv0, distance between ion and point of collision
//				 *tgatomi, SC index of target atom
void	findbc_mjvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             const struct double3D *tgatomlist, const int *Ntgatom,
             const struct ionStruct *ion)
{
	int i;
   double pp, dR;
   struct double3D rIonTg;

   *dRv0=1e20;
   *tgatomi=-1;

   //printf("\nN_SC = %d\n",*Ntgatom);

    // Loop all atoms in tgatomlist to find collision partner
   for (i = 0; (i < *Ntgatom) ; ++i) {

   	//printf("x = %e, y = %e, z = %e\n",
      //  tgatomlist[i].x, tgatomlist[i].y, tgatomlist[i].z);

   	// rIonTg = vector going from ion to target atom i
   	rIonTg.x = tgatomlist[i].x - ion->rUC.x;
   	rIonTg.y = tgatomlist[i].y - ion->rUC.y;
   	rIonTg.z = tgatomlist[i].z - ion->rUC.z;

      // dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
      dR = rIonTg.x * (ion->ev0.x) +
           rIonTg.y * (ion->ev0.y) +
           rIonTg.z * (ion->ev0.z);

      // Check that taget atom
      // 1. is in front of ion
      // 2. is the closest up to now
      if ( (dR > DRMINMIN) && (dR < *dRv0) )
      {
      	// calculate squared impact parameter, pp = p^2
         pp = sqr(rIonTg.x - dR*(ion->ev0.x) ) +   // 6 mult + 5 add
              sqr(rIonTg.y - dR*(ion->ev0.y) ) +
              sqr(rIonTg.z - dR*(ion->ev0.z) );

         // 3. impact parameter is smaller than pmax
         if ( pp < ppmax ) {
            *tgatomi = i;
		      *p = sqrt( pp);
            *dRv0 = dR;
            *R2UC = tgatomlist[i];
         }
      }
   } // for (i ..

   if (*tgatomi != -1)
   {
      // calculate new point of collision, RcUC
		RcUC->x = (*dRv0)*(ion->ev0.x) + ion->rUC.x;
   	RcUC->y = (*dRv0)*(ion->ev0.y) + ion->rUC.y;
   	RcUC->z = (*dRv0)*(ion->ev0.z) + ion->rUC.z;
   }

	// TEST, valid tgatomi: -1 -> NatSC-1
   // removed from code after test period
	if ( ((*tgatomi) < -1) || ((*tgatomi) >= NatSC) )
   {
		printf("Fatal ERROR in function: findtgvib3\n");
   	printf("The program has calculated an invalid target atom index\n");
   	printf("tgatomi = %d, NatSC = %d\n",*tgatomi,NatSC);
      exit(1);
   } // ---
} // end of: findbc_mjvib ---------

//	 ---	findbc_sc_mjvib
// Vibration added before find procedure starts.
// in data:  *ion
// out data: *R2UC, position of closest target atom relative UC
//           *RcUC, point of collision relative UC
//				 *p, imapct parameter
//				 *dRv0, distance between ion and point of collision
//				 *tgatomi, SC index of target atom
//				 *sclist, list of simultaneous collisions
//				 *Nsimcoll, number of simultaneous collisions
void	findbc_sc_mjvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             struct scStruct sclist[], int *Nsimcoll,
             const struct ionStruct *ion)
{
	int i,SCp1,SCp2, potNsimcoll;
   double pp,dR;
   struct double3D rIonTg;
   double u3, phi, Theta;
   struct double3D R2;

   struct bcStruct {
		int tgatomi;
      double dRv0;
		double p;
      struct double3D R2UC;
   }	*bclist;

	double R;

   union rndR_union {
   	double R;           // 64 bit double
      unsigned int i[2];  // i[0] == lo (32) byte of R, i[1] == hi byte of R
	};

   union rndR_union rndR[4];

   *dRv0=1e20;
   *tgatomi=-1;

   *Nsimcoll = 0;

   // Loop all atoms in SC to find collision partner
	SCp1 = SCpointer[ion->SCi];
   SCp2 = SCpointer[(ion->SCi) + 1];

   // Allocate memory for bclist
   // can never be larger than the number of atoms in SubCell
   bclist = (struct bcStruct *) malloc((SCp2-SCp1)*sizeof(struct bcStruct));
   potNsimcoll = 0;

   for (i = SCp1; (i < SCp2) ; ++i) {

      	// --- Add vibration to target atom i at position

   	// Determine four quasi random numbers (here between 1.0 and (2 - 2^-X))
      // from R and random ion.vib numbers.
      // R = |SCatom.r + ion->UCr0|, RR = |R|^2
      R =  sqrt( sqr(SCatom[i].r.x + ion->UCr0.x) +
           		  sqr(SCatom[i].r.y + ion->UCr0.y) +
                 sqr(SCatom[i].r.z + ion->UCr0.z) + EPS );

	   // Construct a positive double rndR with exponent 0 and with
      // the X lowest mantissa bits of R shifted to the X highest
      // bits of rndR. Rlo should then be shiftged (32-1-11 - X) bits
      // and stored in rndR_hi. 0x3ff in rndR_hi unsures positive, 0 exponent
      // double. 10 lowest bits give a shift of 10.
      #define BITSHIFT 10

   	rndR[0].R = (ion->vib.A) * R;
	   rndR[0].i[1] =  ( 0x3ff00000 | ((rndR[0].i[0]<<BITSHIFT) & 0x000fffff) );

   	rndR[1].R = (ion->vib.k.x) * R;
	   rndR[1].i[1] =  ( 0x3ff00000 | ((rndR[1].i[0]<<BITSHIFT) & 0x000fffff) );

   	rndR[2].R = (ion->vib.k.y) * R;
	   rndR[2].i[1] =  ( 0x3ff00000 | ((rndR[2].i[0]<<BITSHIFT) & 0x000fffff) );

   	rndR[3].R = (ion->vib.k.z) * R;
	   rndR[3].i[1] =  ( 0x3ff00000 | ((rndR[3].i[0]<<BITSHIFT) & 0x000fffff) );


      // Calculate normal distributed vibration amplitede u3 (see Xnorm3D)
      // u3 = u30 * sqrt(-2*log(rnd)) * cos(2*pi*rnd)
   	u3 = atom[SCatom[i].id].u3 *
      	  	sqrt( -2.0 * log( rndR[0].R - 1.0 + EPS ) ) *
   		  	cos( 2.0 * pi * (rndR[1].R - 1.0) );

		// 3D random normal vector, evrand, defined by angles: phi and Theta
      // phi = 2*pi * rnd
		phi = 2.0 * pi * (rndR[2].R - 1.0);
      // Theta = acos( 2*rnd - 1)
		Theta = acos( 2.0 * rndR[3].R  - 3.0 );

      // Position of atom including vibration (relative UC)
   	R2.x = SCatom[i].r.x + u3 * sin(Theta)*cos(phi);
   	R2.y = SCatom[i].r.y + u3 * sin(Theta)*sin(phi);
   	R2.z = SCatom[i].r.z + u3 * cos(Theta);
      // ---

   	// rIonTg = vector going from ion to target atom i
   	rIonTg.x = R2.x - ion->rUC.x;
   	rIonTg.y = R2.y - ion->rUC.y;
   	rIonTg.z = R2.z - ion->rUC.z;

      // dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
      dR = rIonTg.x * (ion->ev0.x) +
           rIonTg.y * (ion->ev0.y) +
           rIonTg.z * (ion->ev0.z);


      // Check that taget atom
      // 1. is in front of ion
      //    and a first selection for potential sim. coll.
      if ( (dR > dRmin) && ((dR - *dRv0) < dRsc) )
      {
      	// calculate squared impact parameter, pp = p^2
         pp = sqr(rIonTg.x - dR*(ion->ev0.x) ) +   // 6 mult + 5 add
              sqr(rIonTg.y - dR*(ion->ev0.y) ) +
              sqr(rIonTg.z - dR*(ion->ev0.z) );

         // 2. impact parameter is smaller than pmax
         if ( pp < ppmax ) {
				// Save data for potential sim coll atom in bclist
				bclist[potNsimcoll].R2UC = R2;
				bclist[potNsimcoll].p = sqrt(pp);
				bclist[potNsimcoll].dRv0 = dR;
				bclist[potNsimcoll].tgatomi = i;

            ++potNsimcoll;

            // Update if atom is closest up to now
         	if (dR < *dRv0) {
         		*tgatomi = i;
            	*dRv0 = dR;
               *R2UC = R2;
               *p = sqrt(pp);
            }
         }

      }
   } // for (i ..

   if (*tgatomi != -1)
   {
      // --- Find simultaneous collisions and save sim. coll. data in sclist
		for (i = 0; (i < potNsimcoll); ++i)
      {
      	if (  ((bclist[i].dRv0 - *dRv0) < dRsc ) &&
         		(bclist[i].tgatomi != *tgatomi) )
			// save sim. coll. data in sclist and add thermal vib.
         {
            if (*Nsimcoll > MAXNSIMCOLL){
            	printf("Too many atoms in simultaneous collision!\n Decrease pmax or dRsc.\n");
               exit(1);
            }

            sclist[*Nsimcoll].tgatomi = bclist[i].tgatomi;
				sclist[*Nsimcoll].R2UC = bclist[i].R2UC;
				sclist[*Nsimcoll].p = bclist[i].p;
				sclist[*Nsimcoll].dRv0 = bclist[i].dRv0;


            // Calculate point of collision, RcUC
				sclist[*Nsimcoll].RcUC.x =
            	(sclist[*Nsimcoll].dRv0)*(ion->ev0.x) + ion->rUC.x;
   			sclist[*Nsimcoll].RcUC.y =
            	(sclist[*Nsimcoll].dRv0)*(ion->ev0.y) + ion->rUC.y;
   			sclist[*Nsimcoll].RcUC.z =
            	(sclist[*Nsimcoll].dRv0)*(ion->ev0.z) + ion->rUC.z;

         	++(*Nsimcoll);

         } // end of: if sim. coll.

      } // end of: Find simultaneous collisions

      // Calculate point of collision, *RcUC for closest atom
      RcUC->x = (*dRv0)*(ion->ev0.x) + ion->rUC.x;
      RcUC->y = (*dRv0)*(ion->ev0.y) + ion->rUC.y;
      RcUC->z = (*dRv0)*(ion->ev0.z) + ion->rUC.z;

   }

	// TEST, valid tgatomi: -1 -> NatSC-1
   // removed from code after test period
	if ( ((*tgatomi) < -1) || ((*tgatomi) >= NatSC) )
   {
		printf("Fatal ERROR in function: findtgvib3\n");
   	printf("The program has calculated an invalid target atom index\n");
   	printf("tgatomi = %d, NatSC = %d\n",*tgatomi,NatSC);
      exit(1);
   } // ---

   // Free allocated memory
   free(bclist);

} // end of: findbc_sc_mjvib ---------


// Find target atom randomly distributet in search cylinder.
// in data: *ion
//				 *n, atomic density of target
// out data: *R2UC, position of target atom relative UC
//           *RcUC, point of collision relative UC
//				 *p, imapct parameter
//				 *dRv0, distance between ion and point of collision
//				 *tgatomi, SC index of target atom
void	findrandomtg(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             const struct ionStruct *ion, double *n)
{
   double phi, psinphi, pcosphi;
   struct double3D a3D, e1, e2;


   *p = sqrt(randmj()) * pmax;
   *dRv0 = randmj() * 2.0/( pi * sqr(pmax) * (*n)*1e-24);

   phi = randmj() * 2.0 *pi;

   // Target atom index, tgatomi, is randomly selected by all
   //  atoms in sub-cell strucrure. Change later to atom in Unit Cell
   *tgatomi = (int) floor( randmj()* NatSC );

   // There is a small possibility that randmj() == 1.0
   if (*tgatomi == NatSC) --(*tgatomi);


	RcUC->x = (*dRv0)*(ion->ev0.x) + ion->rUC.x;
   RcUC->y = (*dRv0)*(ion->ev0.y) + ion->rUC.y;
   RcUC->z = (*dRv0)*(ion->ev0.z) + ion->rUC.z;

   // Calculate ortonormal vectors to ion->ev0
   //  a3D = vector that very unlikely will be parallel to ion.ev0
   a3D.x = 17; a3D.y = 2.7183; a3D.z = 0.541667;
   e1 = cross3D(&a3D, &ion->ev0);
   e1 = divk3D(&e1, norm3D(&e1));
   e2 = cross3D(&ion->ev0, &e1);

   psinphi = (*p) * sin(phi);
   pcosphi = (*p) * cos(phi);

   R2UC->x = RcUC->x + e1.x*pcosphi + e2.x*psinphi;
   R2UC->y = RcUC->y + e1.y*pcosphi + e2.y*psinphi;
   R2UC->z = RcUC->z + e1.z*pcosphi + e2.z*psinphi;
} // end of: findrandom target ---------


//	 ---	findbc_robvib
// Vibration added after find procedure completed.
// in data:  *ion
// out data: *R2UC, position of target atom relative UC
//           *RcUC, point of collision relative UC
//				 *p, imapct parameter
//				 *dRv0, distance between ion and point of collision
//				 *tgatomi, SC index of target atom
void	findbc_robvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             const struct ionStruct *ion)
{
	int i,SCp1,SCp2;
   double pp,dR;
   struct double3D x3D, rIonTg;

   *dRv0=1e20;
   *tgatomi=-1;


   // Loop all atoms in SC to find collision partner
	SCp1 = SCpointer[ion->SCi];
   SCp2 = SCpointer[(ion->SCi) + 1];
   for (i = SCp1; (i < SCp2) ; ++i) {

   	// rIonTg = vector going from ion to target atom i
   	rIonTg.x = SCatom[i].r.x - ion->rUC.x;
   	rIonTg.y = SCatom[i].r.y - ion->rUC.y;
   	rIonTg.z = SCatom[i].r.z - ion->rUC.z;

      // dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
      dR = rIonTg.x * (ion->ev0.x) +
           rIonTg.y * (ion->ev0.y) +
           rIonTg.z * (ion->ev0.z);

      // Check that taget atom
      // 1. is in front of ion  more than dRmin
      // 2. is the closest up to now
      if ( (dR > dRmin) && (dR < *dRv0) )
      {
      	// calculate squared impact parameter, pp = p^2
         pp = sqr(rIonTg.x - dR*(ion->ev0.x) ) +   // 6 mult + 5 add
              sqr(rIonTg.y - dR*(ion->ev0.y) ) +
              sqr(rIonTg.z - dR*(ion->ev0.z) );

         // 3. impact parameter is smaller than pmax
         if ( pp < ppmax ) {
            *tgatomi = i;
            *dRv0 = dR;
         }
      }
   } // for (i ..

   if (*tgatomi != -1)
   {
		// Add vibration to target atom *R2UC
      // create random 3D vector, x3D, with normal distributed length
   	Xnorm3D(&x3D);

   	R2UC->x = (SCatom[*tgatomi].r.x + x3D.x * atom[SCatom[*tgatomi].id].u3);
   	R2UC->y = (SCatom[*tgatomi].r.y + x3D.y * atom[SCatom[*tgatomi].id].u3);
   	R2UC->z = (SCatom[*tgatomi].r.z + x3D.z * atom[SCatom[*tgatomi].id].u3);

      // Calculate new dR, P, RcUC
   	// rIonTg = vector going from ion to target atom i
   	rIonTg.x = R2UC->x - ion->rUC.x;
   	rIonTg.y = R2UC->y - ion->rUC.y;
   	rIonTg.z = R2UC->z - ion->rUC.z;

   	// dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
   	*dRv0 = rIonTg.x * (ion->ev0.x) +
      	     rIonTg.y * (ion->ev0.y) +
           	  rIonTg.z * (ion->ev0.z);
      if (*dRv0 < 0.0) *dRv0 = 1e-20;

      // calculate impact parameter p
      *p = sqrt( sqr(rIonTg.x - (*dRv0)*(ion->ev0.x) ) +   // 6 mult + 5 add
                 sqr(rIonTg.y - (*dRv0)*(ion->ev0.y) ) +
                 sqr(rIonTg.z - (*dRv0)*(ion->ev0.z) ) );
      // calculate new point of collision, RcUC
		RcUC->x = (*dRv0)*(ion->ev0.x) + ion->rUC.x;
   	RcUC->y = (*dRv0)*(ion->ev0.y) + ion->rUC.y;
   	RcUC->z = (*dRv0)*(ion->ev0.z) + ion->rUC.z;
   }

	// TEST, valid tgatomi: -1 -> NatSC-1
   // removed from code after test period
	if ( ((*tgatomi) < -1) || ((*tgatomi) >= NatSC) )
   {
		printf("Fatal ERROR in function: findtgvib3\n");
   	printf("The program has calculated an invalid target atom index\n");
   	printf("tgatomi = %d, NatSC = %d\n",*tgatomi,NatSC);
      exit(1);
   } // ---
} // end of: findbc_robvib ---------


//	 ---	findbc_sc_robvib
// Vibration added after find procedure completed.
// in data:  *ion
// out data: *R2UC, position of closest target atom relative UC
//           *RcUC, point of collision relative UC
//				 *p, imapct parameter
//				 *dRv0, distance between ion and point of collision
//				 *tgatomi, SC index of target atom
//				 *sclist, list of simultaneous collisions
//				 *Nsimcoll, number of simultaneous collisions
void	findbc_sc_robvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             struct scStruct sclist[], int *Nsimcoll,
             const struct ionStruct *ion)
{
	int i,SCp1,SCp2, potNsimcoll;
   double pp,dR;
   struct double3D x3D, rIonTg;
   struct bcStruct {
		//struct double3D rIonTg;  added when vibr. is added from the start
		int tgatomi;
      double dRv0;
   }	*bclist;

   *dRv0=1e20;
   *tgatomi=-1;

   *Nsimcoll = 0;

   // --- Loop all atoms in SC to find first collision partner (*dRv0, *tgatomi)
	SCp1 = SCpointer[ion->SCi];
   SCp2 = SCpointer[(ion->SCi) + 1];

   // Allocate memory for bclist
   // can never be larger than the number of atoms in SubCell
   bclist = (struct bcStruct *) malloc((SCp2-SCp1)*sizeof(struct bcStruct));
   potNsimcoll = 0;

   for (i = SCp1; (i < SCp2) ; ++i) {

   	// rIonTg = vector going from ion to target atom i
   	rIonTg.x = SCatom[i].r.x - ion->rUC.x;
   	rIonTg.y = SCatom[i].r.y - ion->rUC.y;
   	rIonTg.z = SCatom[i].r.z - ion->rUC.z;

      // dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
      dR = rIonTg.x * (ion->ev0.x) +
           rIonTg.y * (ion->ev0.y) +
           rIonTg.z * (ion->ev0.z);

      // Check that taget atom
      // 1. is in front of ion  more than dRmin
      //    and a first selection for potential sim. coll.
      if ((dR > dRmin) && ( (dR - *dRv0) < dRsc) )
      {
      	// calculate squared impact parameter, pp = p^2
         pp = sqr(rIonTg.x - dR*(ion->ev0.x) ) +   // 6 mult + 5 add
              sqr(rIonTg.y - dR*(ion->ev0.y) ) +
              sqr(rIonTg.z - dR*(ion->ev0.z) );

         // 2. impact parameter is smaller than pmax
         if ( pp < ppmax ) {
				// Save data for potential sim coll atom in bclist
				// 	bclist[potNsimcoll].rIonTg = rIonTg and  xxxx.pp
            // 	added when vib from start
				bclist[potNsimcoll].dRv0 = dR;
				bclist[potNsimcoll].tgatomi = i;

            ++potNsimcoll;

            // Update if atom is closest up to now
         	if (dR < *dRv0) {
         		*tgatomi = i;
            	*dRv0 = dR;
            }
         }

      }
   } // for (i ..

   if (*tgatomi != -1)
   {
      // --- Find simultaneous collisions and save sim. coll. data in sclist
		for (i = 0; (i < potNsimcoll); ++i)
      {
      	if (  ((bclist[i].dRv0 - *dRv0) < dRsc ) &&
         		(bclist[i].tgatomi != *tgatomi) )
			// save sim. coll. data in sclist and add thermal vib.
         {
            if (*Nsimcoll > MAXNSIMCOLL){
            	printf("Too many atoms in simultaneous collision!\n Decrease pmax or dRsc.\n");
               exit(1);
            }

            sclist[*Nsimcoll].tgatomi = bclist[i].tgatomi;

				// --- Add vibration to target atom *R2UC
		      // create random 3D vector, x3D, with normal distributed length
		   	Xnorm3D(&x3D);

		   	sclist[*Nsimcoll].R2UC.x =
            	(SCatom[bclist[i].tgatomi].r.x +
               x3D.x * atom[SCatom[bclist[i].tgatomi].id].u3);
		   	sclist[*Nsimcoll].R2UC.y =
            	(SCatom[bclist[i].tgatomi].r.y +
               x3D.y * atom[SCatom[bclist[i].tgatomi].id].u3);
		   	sclist[*Nsimcoll].R2UC.z =
            	(SCatom[bclist[i].tgatomi].r.z +
               x3D.z * atom[SCatom[bclist[i].tgatomi].id].u3);

		      // Calculate new dRv0, p, RcUC
		   	// rIonTg = vector going from ion to target atom i
		   	rIonTg.x = sclist[*Nsimcoll].R2UC.x - ion->rUC.x;
		   	rIonTg.y = sclist[*Nsimcoll].R2UC.y - ion->rUC.y;
		   	rIonTg.z = sclist[*Nsimcoll].R2UC.z - ion->rUC.z;

		   	// dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
            dR =	rIonTg.x * (ion->ev0.x) +
                  rIonTg.y * (ion->ev0.y) +
                  rIonTg.z * (ion->ev0.z);
				if (dR < 0.0) dR = 1e-20;

		   	sclist[*Nsimcoll].dRv0 = dR;
	        	// calculate impact parameter p
           	sclist[*Nsimcoll].p =
           		sqrt( sqr(rIonTg.x - dR * (ion->ev0.x) ) +
                 		sqr(rIonTg.y - dR * (ion->ev0.y) ) +
                 		sqr(rIonTg.z - dR * (ion->ev0.z) ) );

		      // calculate new point of collision, RcUC
				sclist[*Nsimcoll].RcUC.x = dR * (ion->ev0.x) + ion->rUC.x;
		   	sclist[*Nsimcoll].RcUC.y = dR * (ion->ev0.y) + ion->rUC.y;
		   	sclist[*Nsimcoll].RcUC.z = dR * (ion->ev0.z) + ion->rUC.z;

         	++(*Nsimcoll);
         } // end of: if sim. coll.

      } // end of: Find simultaneous collisions


		// --- Add vibration to target atom *R2UC
      // create random 3D vector, x3D, with normal distributed length
   	Xnorm3D(&x3D);

   	R2UC->x = (SCatom[*tgatomi].r.x + x3D.x * atom[SCatom[*tgatomi].id].u3);
   	R2UC->y = (SCatom[*tgatomi].r.y + x3D.y * atom[SCatom[*tgatomi].id].u3);
   	R2UC->z = (SCatom[*tgatomi].r.z + x3D.z * atom[SCatom[*tgatomi].id].u3);

      // Calculate new dR, P, RcUC
   	// rIonTg = vector going from ion to target atom i
   	rIonTg.x = R2UC->x - ion->rUC.x;
   	rIonTg.y = R2UC->y - ion->rUC.y;
   	rIonTg.z = R2UC->z - ion->rUC.z;

   	// dR = projection of rIonTg on velocity axis = (rIonTg)dot(ev0)
   	*dRv0 = rIonTg.x * (ion->ev0.x) +
      	     rIonTg.y * (ion->ev0.y) +
           	  rIonTg.z * (ion->ev0.z);
      if (*dRv0 < 0.0) *dRv0 = 1e-20;

      // calculate impact parameter p
      *p = sqrt( sqr(rIonTg.x - (*dRv0)*(ion->ev0.x) ) +   // 6 mult + 5 add
                 sqr(rIonTg.y - (*dRv0)*(ion->ev0.y) ) +
                 sqr(rIonTg.z - (*dRv0)*(ion->ev0.z) ) );
      // calculate new point of collision, RcUC
		RcUC->x = (*dRv0)*(ion->ev0.x) + ion->rUC.x;
   	RcUC->y = (*dRv0)*(ion->ev0.y) + ion->rUC.y;
   	RcUC->z = (*dRv0)*(ion->ev0.z) + ion->rUC.z;
   }

	// TEST, valid tgatomi: -1 -> NatSC-1
   // removed from code after test period
	if ( ((*tgatomi) < -1) || ((*tgatomi) >= NatSC) )
   {
		printf("Fatal ERROR in function: findtgvib3\n");
   	printf("The program has calculated an invalid target atom index\n");
   	printf("tgatomi = %d, NatSC = %d\n",*tgatomi,NatSC);
      exit(1);
   } // ---

   // Free allocated memory
   free(bclist);

} // end of: findbc_sc_robvib ---------




void bcScatter(struct double3D *ev1, double *T, struct double3D *ev2,
               const  struct ionStruct *ion, const struct double3D *R2,
               double M2, const struct double3D *Rcol, double THETA)
{
	double theta, PHI, x;
   struct double3D R2Rc, eR2Rc, e2, a3D;

	// T = E0*4*M1*M2/(M1+M2)^2*sin(THETA/2)^2
	*T =(ion->E)*4.0*(atom[ion->id].M)*M2
   	 /sqr( (atom[ion->id].M) + M2 )*sqr(sin(THETA/2.0));

   // determine laboratory scattering angles from CM-THETA
   x = (atom[ion->id].M)/M2+cos(THETA);

	// The relative precision in the computer is: 2.2204e-16
   if ( fabs(x) < EPS )  //  for example: if THETA == pi and M1 == M2
   	theta = pi/2.0;
   else
	   theta=atan( sin(THETA)/ x );

	PHI=0.5*(pi-THETA);

   // Calculate normalized vector e2 perpendicular
   // to ev0 in the [ev0 x(R2-Rcol)] plane. R2Rc = (R2-Rcol).
   R2Rc = sub3D(R2, Rcol);
	eR2Rc = divk3D(&R2Rc, norm3D(&R2Rc) );

	// e2=(ev0 dot eR2Rc)*ev0-eR2Rc;
   e2 = timesk3D(&ion->ev0, dot3D(&ion->ev0, &eR2Rc) );
   subeq3D(&e2, &eR2Rc );

   // ev1=cos(theta)*ev0+sin(theta)*e2;
   *ev1 = timesk3D(&ion->ev0, cos(theta) );
   a3D = timesk3D(&e2, sin(theta) );
   addeq3D(ev1, &a3D);

	//ev2=cos(PHI)*ev0-sin(PHI)*e2;
   *ev2 = timesk3D(&ion->ev0, cos(PHI) );
   a3D = timesk3D(&e2, sin(PHI) );
   subeq3D(ev2, &a3D);

} // end of: bcScatter ----------------


// Rare event algorithm
void rareEvent(struct ionStruct *ion)
{
	int i;
   struct vibStruct vib0;

   // store original vib data
   vib0 = ion->vib;

   ion->w = ion->w / rareEventMult;

   // Copy (rareEventMult - 1) ions into ion linked list
   for (i = 1; i < rareEventMult; ++i) {

 		// Generate new vib randomnumbers for mj_vib model
   	if (vibalg == MJVIB)
   		mjvibrndnumbers( (&ion->vib));

 		add2ill(ion);
   }

   // restore original vibration for ion(1)
   ion->vib = vib0;
}


      // Add recoiling ion to ion linked list ill


// random number generators  ---------------------------------

void mjvibrndnumbers(struct vibStruct *vib)
// The function creates four random numbers between 1 and RNDAMPLITUDE
// to be used in the mj_vib algorithm
#define RNDAMPLITUDE 10.0
{
		vib->A   = randmj() * (RNDAMPLITUDE - 1.0) + 1.0;
		vib->k.x = randmj() * (RNDAMPLITUDE - 1.0) + 1.0;
		vib->k.y = randmj() * (RNDAMPLITUDE - 1.0) + 1.0;
		vib->k.z = randmj() * (RNDAMPLITUDE - 1.0) + 1.0;
}


void Xnorm3D(struct double3D *x3D)
{

	double phi, Theta, x;

   //   The random normal distributed x is calculated
   //   using the Box-Müller method,
   //   see for example beta math. handbook p. 388
   x = sqrt( -2.0 * log(randmj() + EPS) )
   	 * cos( 2.0 * pi * randmj() );

	// 3D random normal vector, evrand, defined by angles: phi and Theta
	phi = randmj() * 2.0 * pi;
	Theta = acos( 2.0 * randmj() - 1.0 );

	// x3D = x*evrand
   x3D->x = x * sin(Theta)*cos(phi);
   x3D->y = x * sin(Theta)*sin(phi);
   x3D->z = x * cos(Theta);
}  // end of: Xnorm3D

double randmj(void)
// From Mats Hjelm, simplyfied and translated from C++.
// This is a "Minimal Standard" random number generator, adopted for C++
// Descripbed by  Edgar H. Sibley in Communications of the ACM,
// Volume 31 Number 10, October 1988
{
   #define m  2147483647L
   #define q  127773L
   #define a  16807
   #define r  2836

   randomz = a * (long) (randomz % q) - r * (long) (randomz / q);
   if (randomz <= 0) randomz += m;
   return (double) randomz / m;
}


void newrandomseed(void)
{
	 //randomize(); // aparntly only worked for the Borland c-compiler
   srand( (unsigned int) time(0));
   randomz = rand();
}

// ----------------------


double energytest(double Erest)
{
	int i;
   double Etot = 0.0;

   for (i=0; i<Ndistz; i++)
   {
   	Etot += Qa[i];
      Etot += Ta[i];
   };

   return (Etot*dz + Erest) / dose;
}

