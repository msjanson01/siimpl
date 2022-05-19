// file: subcells.c
// Martin Janson, 2003-08-05

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "subcells.h"
#include "mtrxfun.h"
#include "mjmath.h"

#include "mcstructdef.h"
#include "mcglobalextern.h"

#define pi 3.141592654


void cr_subcells( void)
{
	struct scllStruct {
   	struct SCstruct SC;
      struct scllStruct *p;
   };

   struct scllStruct *scll, *tempscll;

   const struct double3D origo = {0.0, 0.0, 0.0};
   struct double3D R0sph[8], R0cyl[12], Rcent[12], R0par[7];
   struct base3D par[7];

	double rmax;
   int NatMC, i, j, ix, iy, iz, count, SCi, inside;
   struct double3D temp3D, ea12, ea23, ea31, a12, a23, a31, r0SC;
   struct int3D nmega;
   struct SCstruct *MCatom;


   // rmax here always the same as pmax
	rmax = pmax;

	// The unitcell is devided into
	// nSC.x * nSC.y * nSC.z subcells determined by the Unit Cell dimensions
   // and the SCsize parametrer
	nSC.x = (int) ceil( norm3D(&(UCbase.a1))/ SCsize);
	nSC.y = (int) ceil( norm3D(&(UCbase.a2))/ SCsize);
	nSC.z = (int) ceil( norm3D(&(UCbase.a3))/ SCsize);

	NSC = nSC.x * nSC.y * nSC.z;

	// Calculate subcell base and transformation base
   SCbase.a1 = divk3D(&UCbase.a1, nSC.x);
   SCbase.a2 = divk3D(&UCbase.a2, nSC.y);
   SCbase.a3 = divk3D(&UCbase.a3, nSC.z);

   on2SCbase = invbase(SCbase);

 	// ----- 	Allocate memory for Global variable SCpointer ----------------
   SCpointer = (int *) malloc((NSC + 1)*sizeof(int));


	// Calculate subcell base vectors
	SCbase.a1  = divk3D( &(UCbase.a1), nSC.x);
	SCbase.a2  = divk3D( &(UCbase.a2), nSC.y);
	SCbase.a3  = divk3D( &(UCbase.a3), nSC.z);


	// Calculate normal vectors to the unitcell (and SC) surfaces

  temp3D = cross3D( &(UCbase.a1), &(UCbase.a2));
	ea12 = divk3D( &temp3D, norm3D( &temp3D) );

  temp3D = cross3D( &(UCbase.a2), &(UCbase.a3));
  ea23 = divk3D( &temp3D, norm3D( &temp3D) );

  temp3D = cross3D( &(UCbase.a3), &(UCbase.a1));
	ea31 = divk3D( &temp3D, norm3D( &temp3D) );


	// Create megacell (MC) so that all points within
	// a distance rmax from the unitcell surface
	// are hosted within MC.
	// MC is constructed by (2*n1mega+1)*(2*n2mega+1)*(2*n3mega+1) unitcells.
   // the atom data (SCstruct) are stored in MCatom;

	nmega.x = (int) ceil(rmax/dot3D(&(UCbase.a1), &ea23));
	nmega.y = (int) ceil(rmax/dot3D(&(UCbase.a2), &ea31));
	nmega.z = (int) ceil(rmax/dot3D(&(UCbase.a3), &ea12));

	NatMC = NatUC*(1+2*nmega.x)*(1+2*nmega.y)*(1+2*nmega.z);

	// 	Allocate memory for MC structure
   MCatom = (struct SCstruct *) malloc( NatMC * sizeof(struct SCstruct));

   count = 0;
   for (iz = -nmega.z; iz <= nmega.z; iz++)
   	for (iy = -nmega.y; iy <= nmega.y; iy++)
   		for (ix = -nmega.x; ix <= nmega.x; ix++)
        for (i=0; i < NatUC; i++) {
          	MCatom[count] = UCatom[i];
            // Calculate reference point of unit cell
        	  temp3D.x = ix;
            temp3D.y = iy;
				    temp3D.z = iz;
            temp3D = vecttransf(temp3D, UCbase);
            // and add it to the position of the UC atom
            addeq3D(&(MCatom[count].r), &temp3D);
            ++count;
         };


	//   Find atoms from the MC that belongs to each subcell, i.e
	// the atoms that are within a distance rmax from the
	// subcell limiting surface. This volume can be represented by
	// 8 spheres (r=rmax), 12 cylinders (r=rmax) and 6+1 parallel epided.

	// Define R0 of the 8 spheres (relative to SC).
   R0sph[0] = origo;
   R0sph[1] = SCbase.a1;
   R0sph[2] = add3D( &(SCbase.a1), &(SCbase.a2) );
   R0sph[3] = SCbase.a2;
   R0sph[4] = SCbase.a3;
   R0sph[5] = add3D( &(SCbase.a1), &(SCbase.a3) );
   // R06 = a1 + a2 + a3
   R0sph[6] = add3D( &(SCbase.a1), &(SCbase.a2) );
   addeq3D( &R0sph[6], &(SCbase.a3));
   R0sph[7] = add3D( &(SCbase.a2), &(SCbase.a3) );

	// Define R0 and Rcent of the 12 cylinders (relative to SC).
	R0cyl[0] = origo;
	R0cyl[1] = SCbase.a1;
	R0cyl[2] = origo;
	R0cyl[3] = SCbase.a2;
	R0cyl[4] = SCbase.a3;
   R0cyl[5] = add3D( &(SCbase.a1), &(SCbase.a3) );
	R0cyl[6] = SCbase.a3;
   R0cyl[7] = add3D( &(SCbase.a2), &(SCbase.a3) );
	R0cyl[8] = origo;
	R0cyl[9] = SCbase.a1;
   R0cyl[10] = add3D( &(SCbase.a1), &(SCbase.a2) );
	R0cyl[11] = SCbase.a2;

   Rcent[0] = SCbase.a1;
   Rcent[1] = SCbase.a2;
   Rcent[2] = SCbase.a2;
   Rcent[3] = SCbase.a1;
   Rcent[4] = SCbase.a1;
   Rcent[5] = SCbase.a2;
   Rcent[6] = SCbase.a2;
   Rcent[7] = SCbase.a1;
   Rcent[8] = SCbase.a3;
   Rcent[9] = SCbase.a3;
   Rcent[10] = SCbase.a3;
   Rcent[11] = SCbase.a3;

   // Define the 6+1 paralell epipeds
   a12 = timesk3D(&ea12, rmax);
   a23 = timesk3D(&ea23, rmax);
   a31 = timesk3D(&ea31, rmax);

 	R0par[0] = timesk3D(&a31, -1);  // -a31
   R0par[1] = SCbase.a1;
   R0par[2] = SCbase.a2;
   R0par[3] = timesk3D(&a23, -1);  // -a23
   R0par[4] = timesk3D(&a12, -1);  // -a12
   R0par[5] = SCbase.a3;
   R0par[6] = origo;

	par[0] = vects2base(SCbase.a1, a31, SCbase.a3);
   par[1] = vects2base(SCbase.a2, SCbase.a3, a23);
   par[2] = vects2base(SCbase.a1, a31, SCbase.a3);
   par[3] = vects2base(SCbase.a2, SCbase.a3, a23);
   par[4] = vects2base(SCbase.a1, SCbase.a2, a12);
   par[5] = vects2base(SCbase.a1, SCbase.a2, a12);
   par[6] = SCbase;

   // --------------------------------

	// Initiate SC linked list (scll)
   scll = (struct scllStruct *) malloc(sizeof(struct scllStruct));
   // bottom of stack is marked with ill->ion.id = -1
   scll->SC.id = -1;

   // Set the first SCpointer to 0
   SCpointer[0] = 0;

   NatSC = 0;
	SCi = 0; // counter of sub cell
 	for (iz=1; iz <= nSC.z; iz++)
		for (iy=1; iy <= nSC.y; iy++)
			for (ix=1; ix <= nSC.x; ix++) {

     			// base vector r0SC for sub cell SCi
            temp3D.x = ix - 1.0;
            temp3D.y = iy - 1.0;
            temp3D.z = iz - 1.0;
      		r0SC = vecttransf( temp3D, SCbase);

				// Find atoms belonging to sub cell SCi from the MC atoms list
		      for (i=0; i < NatMC; i++) {
         		inside = 0;

         		// look in 8 spheres
         		j = 0;
               while ( (j < 8) & (!inside) ) {
               	temp3D = add3D(&R0sph[j], &r0SC);
               	inside = insideSphere(MCatom[i].r, temp3D, rmax);
            		++j;
         		};

         		// look in 12 cylinders
         		j = 0;
               while ( (j < 12) & (!inside) ) {
               	temp3D = add3D(&R0cyl[j], &r0SC);
               	inside = insideCyl(MCatom[i].r, temp3D, Rcent[j], rmax);
            		++j;
					};

         		// look in 6+1 paralell epideds
         		j = 0;
               while ( (j < 7) & (!inside) ) {
               	temp3D = add3D(&R0par[j], &r0SC);
               	inside = insideParep(MCatom[i].r, temp3D, par[j]);
            		++j;
					};

         		if (inside) {
               	// Add MC atom to sub cell linked list
                  tempscll = scll;
  					   if ( ( scll = (struct scllStruct *)
                  	   malloc(sizeof(struct scllStruct)) )== NULL ) {
					   	printf("Out of memory in add2ill function!\n");
      					exit(1);
   					}
					   scll->p = tempscll;
					   scll->SC = MCatom[i];

                  ++NatSC;
					};
            }; // end of for i: MCatom loop
         // update SCpointer[SCi]
         ++SCi;
         SCpointer[SCi] = NatSC;

      }; // end of for ix, iy, iz


	// Free allocated memory of the MCatom array
   free(MCatom);

 	// Allocate memory for Global variable SCatom ----------------
   SCatom = (struct SCstruct *) malloc(NatSC*sizeof(struct SCstruct));

   for (i = NatSC; i > 0; i--) {
   	SCatom[i - 1] = scll->SC;

      // remove top of the stack and change pointer scll
      tempscll = scll->p;
      free(scll);
      scll = tempscll;
   };

   //Check that scll now points at first cell withh .id = -1
   if (scll->SC.id != -1) {
		printf("Last scll->SC.id != -1\n");
		exit(1);
	}

   // Free last cell of scll list
   free(scll);

}; // end of cr_subcells




// ---  Special functions used by crSubCells

int insideSphere(struct double3D Rp, struct double3D R0,
					  double radius)
{
	struct double3D temp3D;

   temp3D = sub3D(&Rp, &R0);
   if (norm3D( &temp3D ) < radius )
   	return 1;
   else
   	return 0;

}  //end of: insideSphere  ------------


int insideCyl(struct double3D Rp, struct double3D R0,
				  struct double3D Rcent, double radius)
{
   double nRc, r, l;
   struct double3D er, temp3D;

	nRc = norm3D(&Rcent);

	er = divk3D( &Rcent, nRc);

	Rp = sub3D( &Rp, &R0);

   temp3D = cross3D( &er, &Rp);
	r = norm3D(&temp3D);

	l = dot3D( &er, &Rp);

	if ((r <= radius) && (l >= 0) && (l <= nRc))
   	return 1;
   else
   	return 0;

}  // end of: insideCyl --------------------


int insideParep(struct double3D Rp, struct double3D R0,
					 struct base3D a)
{
// a1,a2,a3 (in a) must be defined according to the righ hand rule!
   struct double3D i3D;
  // tempbase3D = invbase(a);

	Rp = sub3D( &Rp, &R0);

   i3D = floorvecttransf( Rp, invbase(a));

   if ( ( (int) i3D.x == 0) && ( (int) i3D.y  == 0) && ( (int) i3D.z  == 0))
   	return 1;
   else
   	return 0;

} // end of: insideParep ------------------           */
