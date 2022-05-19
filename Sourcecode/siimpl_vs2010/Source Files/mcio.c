// file: mcio.c
// Martin Janson, 2000-03-09

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mcio.h"
#include "mcstd.h"
#include "mtrxfun.h"
#include "mjmath.h"
#include "spstd.h"
                                                                         
#include "mcstructdef.h"
#include "mcglobalextern.h"                                      


#define pi 3.141592654





void structureSize(char *filename, int *NatUC,
											  int *Nspecies, int *Ndz)
// finds structure size data in *filename
{
	FILE *fid;
   char firstStr[100];
	float x,y,z;
   int n;
   int i;

   if ((fid=fopen(filename, "r")) == NULL) {
		printf("Can not open: %s\n",filename);
    	exit(1);
   }

    while (!feof(fid))  {
      firstStr[0]='\0';
   	fscanf(fid,"%s",firstStr);

      // read simfile, (recursive function call)
      if (!strcmp(firstStr,"readsimfile") )
      {
      	fscanf(fid,"%s",firstStr);
         structureSize(firstStr, NatUC, Nspecies, Ndz);
      }

      // Number of atoms in UC---------------------------
      else if (!strcmp(firstStr,"UCatoms")) {
      	fskipline(fid);   // Reads away line with "UCatoms"

			i=0;
			while (fscanf(fid,"%d%e%e%e",&n,&x,&y,&z) == 4) {

            fskipline(fid);
            ++i;
         }

			*NatUC = i;
         fskipline(fid);
      } // 	--------------------------------------------


      // Numbber of target atom species
     	else if (!strcmp(firstStr,"Z2")) {
         i=0;
         while (fscanf(fid,"%d",&n)==1) {
            ++i;
         };
			*Nspecies = i;
		}

      // Number of slabs in output distributions
   	else if (!strcmp(firstStr,"Ndistz")) {
      	fscanf(fid,"%d",Ndz);
     	}
      // read away rest of the line until CR OR(!) LF --------------------
      fskipline(fid);
   }

   fclose(fid);
}  // structureSize ----------------------------------------------------



int readSimData(char *filename)
// Reads simulation data from file "filename",
{
	FILE *fid;
   char firstStr[100],secStr[100];
	float x,y,z;
   struct double3D dbl3Dtmp;
   int n,i,j;

   if ((fid=fopen(filename, "r")) == NULL) {
		printf("Can not open: %s\n",filename);
    	exit(1);
   }

    while ( (fscanf(fid,"%s",firstStr)==1) ) {



      // read simfile, (recursive function call)  ----
      if (!strcmp(firstStr,"readsimfile") ){
      	if (fscanf(fid,"%s",firstStr)==1) readSimData(firstStr);
         else exitRD("readsimfile");
      }

      // Import Czion, CzI, and CzV files
      else if (!strcmp(firstStr,"importCzion") ){
      	if (fscanf(fid,"%s",secStr)!=1) exitRD("importCzion");

         if (!strcmp(secStr,"current") )
         	secStr[0] = 0;

			importCzion(secStr);
      }
      else if (!strcmp(firstStr,"importCzI") ){
      	if (fscanf(fid,"%s",secStr)!=1) exitRD("importCzI");

         if (!strcmp(secStr,"current") )
         	secStr[0] = 0;

			importCzI(secStr);
      }
      else if (!strcmp(firstStr,"importCzV") ){
      	if (fscanf(fid,"%s",secStr)!=1) exitRD("importCzV");

         if (!strcmp(secStr,"current") )
         	secStr[0] = 0;

			importCzV(secStr);
      }

		// Lattice structure -------------------------------

      //  Unit Cell
   	else if (!strcmp(firstStr,"UCbase")) {
      	n = 0;

      	fskipline(fid);
      	n += fscanf(fid,"%e%e%e",&x,&y,&z);
			UCbase.a1.x=x; UCbase.a1.y=y; UCbase.a1.z=z;

      	fskipline(fid);
      	n += fscanf(fid,"%e%e%e",&x,&y,&z);
			UCbase.a2.x=x; UCbase.a2.y=y; UCbase.a2.z=z;

      	fskipline(fid);
      	n += fscanf(fid,"%e%e%e",&x,&y,&z);
			UCbase.a3.x=x; UCbase.a3.y=y; UCbase.a3.z=z;

         if (n!= 9)  exitRD("UCbase");
      }
     	else if (!strcmp(firstStr,"Z2")) {
         i=1;
         while (fscanf(fid,"%d",&n)==1) {
         	atom[i].Z=n;
            ++i;
         };
		}
     	else if (!strcmp(firstStr,"M2")) {
         i=1;
         while (fscanf(fid,"%e",&x)==1) {
         	atom[i].M=x;
            ++i;
         };
         if (i > Nspecies+1) exitRD2many("M2");
		}

      else if (!strcmp(firstStr,"UCatoms")) {
      	fskipline(fid);   // Reads away line with "UCatoms"

			i=0;
			while (fscanf(fid,"%d%e%e%e",&n,&x,&y,&z) == 4) {

            UCatom[i].id = n,
            UCatom[i].r.x = x,
            UCatom[i].r.y = y,
            UCatom[i].r.z = z;

            fskipline(fid);
            ++i;
         }

			NatUC = i;
         fskipline(fid);
      } // 	--------------------------------------------

      // Target atom data ---------------
   	else if ( !strcmp(firstStr,"randomtarget") ||
               (!strcmp(firstStr,"randomtarget;")) ){
      	randomtarget = 1;
     	}

   	else if (!strcmp(firstStr,"natomamorph")) {
      	if (fscanf(fid,"%e",&x)==1)  nAtomAmorph = x; //cm-3
         else exitRD("natomamorph");
         // atomic density in amorphous surface layer, nAtomAmorph (cm-3)
     	}

     	else if (!strcmp(firstStr,"Edispl")) {
         i=1;
         while (fscanf(fid,"%e",&x)==1) {
         	atom[i].Ed=x;
            ++i;
         };
         if (i > Nspecies+1) exitRD2many("Edispl");
		}
     	else if (!strcmp(firstStr,"vibu1")) {
         i=1;
         while (fscanf(fid,"%e",&x)==1) {
         	atom[i].u3=sqrt(3.0)*x;
            ++i;
        };
        if (i > Nspecies+1) exitRD2many("vibu1");
		}
 		// -----------------------------------------------------------

      // SC size definition
   	else if (!strcmp(firstStr,"SCsize")) {
      	if (fscanf(fid,"%e",&x)==1)  SCsize = x; // Å
         else exitRD("SCsize");
     	}


      // implanted ion data ----------------------------------------
   	else if (!strcmp(firstStr,"Z1")) {
      	if (fscanf(fid,"%d",&n)==1) atom[0].Z=n;
         else exitRD("Z1");
      }
   	else if (!strcmp(firstStr,"M1")) {
      	if (fscanf(fid,"%e",&x)==1) atom[0].M=x;
         else exitRD("M1");
      }
   	else if (!strcmp(firstStr,"E0")) {
      	if (fscanf(fid,"%e",&x)==1) ion0.E=x;
         else exitRD("E0");
      }
   	else if (!strcmp(firstStr,"impdir")) {
      	if (fscanf(fid,"%e%e%e",&x,&y,&z)==3) {
				// ensure normalisation
  				ion0.ev0.x=x/sqrt(x*x+y*y+z*z);
  				ion0.ev0.y=y/sqrt(x*x+y*y+z*z);
  				ion0.ev0.z=z/sqrt(x*x+y*y+z*z);
         }
         else exitRD("impdir");
      }
   	else if (!strcmp(firstStr,"impdirUC")) {
      	if (fscanf(fid,"%e%e%e",&x,&y,&z)==3) {

            dbl3Dtmp.x = x;
            dbl3Dtmp.y = y;
            dbl3Dtmp.z = z;

            dbl3Dtmp = vecttransf( dbl3Dtmp, UCbase);

            x = (float) dbl3Dtmp.x;
            y = (float) dbl3Dtmp.y;
            z = (float) dbl3Dtmp.z;

				// ensure normalisation
  				ion0.ev0.x = x/sqrt(x*x+y*y+z*z);
  				ion0.ev0.y = y/sqrt(x*x+y*y+z*z);
  				ion0.ev0.z = z/sqrt(x*x+y*y+z*z);
         }
         else exitRD("impdirUC");
      }
   	else if (!strcmp(firstStr,"impdirrotyrotz")) {
      // i.e. first rotate around ey, then rotate around ez
      	if (fscanf(fid,"%e%e",&y,&x)==2) {
  				ion0.ev0.x = cos(x*pi/180) * sin(y*pi/180);
         	ion0.ev0.y = sin(x*pi/180) * sin(y*pi/180);
         	ion0.ev0.z = cos(y*pi/180);
         }
         else exitRD("impdirrotyrotz");
      }
   	else if (!strcmp(firstStr,"dose")) {
      	if (fscanf(fid,"%e",&x)==1) dose=x;
         else exitRD("dose");
      }


      // definition of sample surface -----------------
   	else if (!strcmp(firstStr,"Rsurf")) {
      	if (fscanf(fid,"%e%e%e",&x,&y,&z)==3) {
	  			Rsurf.x=x;
            Rsurf.y=y;
            Rsurf.z=z;
         }
         else exitRD("Rsurf");
      }
   	else if (!strcmp(firstStr,"surfnorm")) {
      	if (fscanf(fid,"%e%e%e",&x,&y,&z)==3) {
	  			// ensure normalisation
  				ensurf.x=x/sqrt(x*x+y*y+z*z);
  				ensurf.y=y/sqrt(x*x+y*y+z*z);
  				ensurf.z=z/sqrt(x*x+y*y+z*z);
         }
         else exitRD("surfnorm");
      }
   	else if (!strcmp(firstStr,"surfnormUC")) {
      	if (fscanf(fid,"%e%e%e",&x,&y,&z)==3) {

            dbl3Dtmp.x = x;
            dbl3Dtmp.y = y;
            dbl3Dtmp.z = z;

            dbl3Dtmp = vecttransf( dbl3Dtmp, UCbase);

            x = (float) dbl3Dtmp.x;
            y = (float) dbl3Dtmp.y;
            z = (float) dbl3Dtmp.z;

				// ensure normalisation
  				ensurf.x = x/sqrt(x*x+y*y+z*z);
  				ensurf.y = y/sqrt(x*x+y*y+z*z);
  				ensurf.z = z/sqrt(x*x+y*y+z*z);
         }
         else exitRD("surfnormUC");
      }
   	else if (!strcmp(firstStr,"surfnormrotyrotz")) {
      // i.e. first rotate around ey, then rotate around ez
      	if (fscanf(fid,"%e%e",&y,&x)==2) {
  				ensurf.x = cos(x*pi/180) * sin(y*pi/180);
         	ensurf.y = sin(x*pi/180) * sin(y*pi/180);
         	ensurf.z = cos(y*pi/180);
         }
         else exitRD("surfnormrotyrotz");
      }
   	else if (!strcmp(firstStr,"layerwidth")) {
      	if (fscanf(fid,"%e",&x)==1)  layerWidth = x;
         else exitRD("layerwidth");
      }

     else if (!strcmp(firstStr,"randsurflayer")) {
      	if (fscanf(fid,"%e",&x)==1) randsurflayer = x;
         else exitRD("randsurflayer");
      }


		// definition of implanted area   ----------------
   	else if (!strcmp(firstStr,"Rimplant")) {
      	if (fscanf(fid,"%e%e%e",&x,&y,&z)==3) {
	  			Rimp.x=x; Rimp.y=y; Rimp.z=z;
         }
         else if (fscanf(fid,"%s",&secStr)==1)  {
         	if ( (!strcmp(secStr,"Rsurf"))||(!strcmp(secStr,"Rsurf;")) ){
            	Rimp = Rsurf;
            } else exitRD("Rimplant");
         }
      }
   	else if (!strcmp(firstStr,"implantarea")) {
      	if (fscanf(fid,"%e%e",&x,&y)==2) {
				impArea1 = x;
   	      impArea2 = y;
         }
			else exitRD("implantarea");
      }

   	else if (!strcmp(firstStr,"surfaceenergy")) {
      	if (fscanf(fid,"%e",&x)==1)  U_s = x;
         else exitRD("surfaceenergy");
      }

      // nuclear stopping data
   	else if (!strcmp(firstStr,"IAPfunction")) {
      	if (fscanf(fid,"%d",&n)==1)  IAPfunction = n;
         else exitRD("IAPfunction");
      }

      // electronic stopping power data ---------------
   	else if ( (!strcmp(firstStr,"Efermi")) ||
       			(!strcmp(firstStr,"EFERMI")) ){
      	if (fscanf(fid,"%e",&x)==1)  E_fermi = x;
         else exitRD("Efermi");
      }

     	else if  (!strcmp(firstStr,"esp.c") ){
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
         	if (i > Nspecies)  exitRD2many("esp.c");
         	// Low velocity Se according to LS and Braggs rule
   			esp[i].A1 = 0.0;
      		for (j=1; j<=Nspecies; j++) {
            esp[i].A1 +=
         		(atom[j].x) * SeLS(1.0, atom[i].Z, atom[i].M, atom[j].Z);
      		}
      		esp[i].A1 *=  x;
       		++i;
         };
  		}
     	else if (!strcmp(firstStr,"esp.p"))    {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
  		     	if (i > Nspecies)  exitRD2many("esp.p");
           	esp[i].p = x;
            ++i;
         };
		}
     	else if ( !strcmp(firstStr,"esp.p2" ) ){
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
  		     	if (i > Nspecies)  exitRD2many("esp.p2");
          		esp[i].p2 = x;
         	++i;
         };
		}
     	else if ( !strcmp(firstStr,"esp.A1" )) {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
         	if (i > Nspecies)  exitRD2many("esp.A1");
        		esp[i].A1 = x * 1.0e-2; //    eV/(10^15cm^-2) -> keV*Å^2
         	++i;
         };
  		}
     	else if ( !strcmp(firstStr,"esp.A2" )) {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
         	if (i > Nspecies)  exitRD2many("esp.A2");
        		esp[i].p = x;
        		esp[i].p2 = 0.0;
         	++i;
         };
  		}
      // NOTE that the (E/1000) in the Se_Hi formula is performed
      // directly on the A3-A5 parameters here
     	else if ( !strcmp(firstStr,"esp.A3" )) {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
         	if (i > Nspecies)  exitRD2many("esp.A3");
        		esp[i].A3 = x * 1.0e1; // MeV*eV/(10^15cm^-2) -> keV^2*Å^2
         	++i;
         };
		}
     	else if ( !strcmp(firstStr,"esp.A4" )) {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
       	if (i > Nspecies)  exitRD2many("esp.A4");
          		esp[i].A4 = x * 1e3;    // MeV^-1 -> keV^-1
         	++i;
         };
		}
     	else if ( !strcmp(firstStr,"esp.A5" )) {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
       	if (i > Nspecies)  exitRD2many("esp.A5");
          		esp[i].A5 = x * 1e-3;   // MeV -> keV
         	++i;
         };
		}
     	else if ( !strcmp(firstStr,"esp.fl" ) ){
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
  		     	if (i > Nspecies)  exitRD2many("esp.fl");
          		esp[i].lf = x;
         	++i;
         };
		}
     	else if (!strcmp(firstStr,"esp.s"))  {
         i=0;
         while (fscanf(fid,"%e",&x)==1) {
  		     	if (i > Nspecies)  exitRD2many("esp.s");
          		esp[i].s = x;
         	++i;
         };
		}

		// MC statistics  ------------------------

   	else if (!strcmp(firstStr,"Nions")) {
      	if (fscanf(fid,"%d",&n)==1)  Nions = n;
         else exitRD("Nions");
      }
   	else if (!strcmp(firstStr,"Ndistz")) {
      	if (fscanf(fid,"%d",&n)==1)  Ndistz = n;
         else exitRD("Ndistz");
      }
   	else if ( (!strcmp(firstStr,"saveionsR3"))
      			|| (!strcmp(firstStr,"saveionsR3;")) ){
      	saveionsR3 = 1;
      }
   	else if ( (!strcmp(firstStr,"saveIVR3"))
  					|| (!strcmp(firstStr,"saveIVR3;")) ){
      	saverecoilsR3 = 1;
      }
     	else if (!strcmp(firstStr,"saveatdose")){
         i=0;
         while ((fscanf(fid,"%e",&x)==1) && (i<10)) {
         	savedose[i] = x;
            ++i;
         };
         Nsavedose = i;
		}
   	else if ( (!strcmp(firstStr,"calcNEP"))
  					|| (!strcmp(firstStr,"calcNEP;")) ){
      	calcFnep = 1;
      }

      // Rare event model
   	else if (!strcmp(firstStr,"rareeventmult")) {
      	if (fscanf(fid,"%e",&x)==1) {
         	rareEventMult = (int) x;
         }
         else exitRD("rareeventmult");
      }

     	else if (!strcmp(firstStr,"raredepth")){
         i=0;
         while ((fscanf(fid,"%e",&x)==1) && (i < MAXNRAREDEPTH)) {
         	rareDepth[i] = x;
            ++i;
         };
         if (i == MAXNRAREDEPTH) {
         	printf("Maximum number of rare depths is %d\n", MAXNRAREDEPTH);
            exit(1);
         };
         NrareDepth = i;

		}



      // BCA parameters  -----------------------------
   	else if (!strcmp(firstStr,"pmax")) {
      	if (fscanf(fid,"%e",&x)==1) {
         	pmax=x;
         }
         else exitRD("pmax");
      }
   	else if (!strcmp(firstStr,"Estop")) {
      	if (fscanf(fid,"%e",&x)==1) Estop=x;
         else exitRD("Estop");
      }
      else if ( (!strcmp(firstStr,"randomize"))
      			|| (!strcmp(firstStr,"randomize;")) ){
      	newrandomseed();
      }
   	else if (!strcmp(firstStr,"dRsc")) {
      	if (fscanf(fid,"%e",&x)==1) {
         	dRsc = x;
         }
         else exitRD("dRsc");
      }
      else if ( (!strcmp(firstStr,"simultcollision"))
      			 || (!strcmp(firstStr,"simultcollision;")) ){
      	simcollalg = 1;
      }
   	else if (!strcmp(firstStr,"vibalg")) {
      	if (fscanf(fid,"%d",&n)==1)  vibalg = n;
         else exitRD("vibalg");
      }

      // damage model   ----------------------------
     	else if  (!strcmp(firstStr,"damage.cascade")) {
      	if (fscanf(fid,"%d",&n)== 1)  recoil = n;
         else exitRD("damage.cascade");
      }
   	else if (!strcmp(firstStr,"damage.c_a")) {
      	if (fscanf(fid,"%e",&x)==1) rndscatfac = x;
         else exitRD("damage.c_a");
      }


      // read comment signs
      else if ( (firstStr[0]=='%') ||
      			 (firstStr[0]=='$')    );

      else {
      	printf(" '%s' is not a recognized instruction!\n",firstStr);
         printf("Check instruction manual for correct format!\n");
		   exit(1);
      }

      // read away rest of the line until CR OR(!) LF --------------------
      fskipline(fid);
   }

   fclose(fid);
	return 1;
}  // end of: readSimData ----------------------------------------------------




// Error message functions
void exitRD(char *mesg)
{
	printf("Erroneous format when specifying '%s' !\n",mesg);
   printf("Check instruction manual for correct format!\n");
   exit(1);
}
void exitRD2many(char *mesg)
{
	printf("Wrong number of input values in '%s' for current SC structure!\n"
   	,mesg);
   printf("Check instruction manual for correct format!\n");
   exit(1);
}
// -------------------------------------------------------------------

//  ---- Import Cz data functions

void 	importCzion(const char path[100])
{
	FILE *fid;
   char filename[115];
	const char name[] = "Czion.txt";
   float x, y;
   int i;

	// Create filename = path + name
   strcpy(filename, path);
   strcpy(&filename[strlen(filename)], name);

   if ((fid=fopen(filename, "r")) == NULL) {
		printf("Can not open: %s (importCzion)\n",filename);
    	exit(1);
   }

   // Read away first header line
   fskipline(fid);

	i = 0;
  	while (fscanf(fid,"%e%e", &x, &y) == 2) {
		Czion[i] = y;
		++i;
	}
   fclose(fid);

   if (i != Ndistz) {
    	printf("The dimension of the imported 'Czion.txt' file is "
      "not compatible with the 'Ndistz' of the current simulation!\n");
       exit(1);
   }
} // end of importCzion function

void 	importCzV(const char path[100])
{
	FILE *fid;
   char filename[115];
	const char name[] = "CzV.txt";
   float x;
   int i, j;

	// Create filename = path + name
   strcpy(filename, path);
   strcpy(&filename[strlen(filename)], name);

   if ((fid=fopen(filename, "r")) == NULL) {
		printf("Can not open: %s (importCzV)\n",filename);
    	exit(1);
   }

   // Read away first header line
   fskipline(fid);

	i = 0;
  	while (fscanf(fid,"%e", &x) == 1) {
		for (j=0; j < Nspecies; j++)
			if (fscanf(fid,"%e", &x) == 1)
         	CzV[i + j*Ndistz] = x;
         else
         	exitRD("importCzV");
		++i;
	}
   fclose(fid);

   if (i != Ndistz) {
    	printf("The dimension of the imported 'CzV.txt' file is "
      "not compatible with the 'Ndistz' of the current simulation!\n");
       exit(1);
   }
} // end of importCzV function



void 	importCzI(const char path[100])
{
	FILE *fid;
   char filename[115];
	const char name[] = "CzI.txt";
   float x;
   int i, j;

	// Create filename = path + name
   strcpy(filename, path);
   strcpy(&filename[strlen(filename)], name);

   if ((fid=fopen(filename, "r")) == NULL) {
		printf("Can not open: %s (importCzI)\n",filename);
    	exit(1);
   }

   // Read away first header line
   fskipline(fid);

	i = 0;
  	while (fscanf(fid,"%e", &x) == 1) {
		for (j=0; j < Nspecies; j++)
			if (fscanf(fid,"%e", &x) == 1)
         	CzI[i + j*Ndistz] = x;
         else
         	exitRD("importCzI");
		++i;
	}
   fclose(fid);

   if (i != Ndistz) {
    	printf("The dimension of the imported 'CzI.txt' file is "
      "not compatible with the 'Ndistz' of the current simulation!\n");
       exit(1);
   }
} // end of importCzI function




//----------     Output file functions     ----------------------------

void saveCzData(int action)
// action = 0; open files
//				1; write to files
//				2; close files
{
	static FILE *fCzion, *fQTz,*fCzI, *fCzV, *fp;
   int i,j;

   switch (action){
	   case 0: // Open files
		   if ((fCzion = fopen("Czion.txt", "w")) == NULL) {
   			printf("Can not open: Czion.txt for writing\n");
            exit(1);
   		};
		   if ((fQTz = fopen("QTz.txt", "w")) == NULL) {
   			printf("Can not open: QTz.txt for writing\n");
            exit(1);
   		};
		   if ((fCzI = fopen("CzI.txt", "w")) == NULL) {
   			printf("Can not open: CzI.txt for writing\n");
            exit(1);
   		};
		   if ((fCzV = fopen("CzV.txt", "w")) == NULL) {
   			printf("Can not open: CzV.txt for writing\n");
            exit(1);
   		};
		   if ((fp = fopen("NEP.txt", "w")) == NULL) {
   			printf("Can not open: NEP.txt for writing\n");
            exit(1);
   		};

	   	fprintf(fCzion, "%% Depth (A)\tConc. (cm-3)\n");
	   	fprintf(fQTz, "%% Depth (A)\tQ (keV/Å/cm2)\tT (keV/A/cm2)\n");
   		fprintf(fCzI, "%% Depth (A)\tI1 (cm^-3)\tI2 (cm^-3)\t...\n");
   		fprintf(fCzV, "%% Depth (A)\tV1 (cm^-3)\tV2 (cm^-3)\t...\n");
   		fprintf(fp, "%% Depth (A)\tNEP1 (Å-1)\tNEP2 (Å-1)\t...\n");

         break;

      case 1: // Write to files
		   for (i = 0; i<Ndistz; i++)
         {

   			fprintf(fCzion, "%e\t%e\n", (dz*i+dz/2.0), Czion[i]);
   			fprintf(fQTz, "%e\t%e	%e\n", (dz*i+dz/2.0), Qa[i], Ta[i]);

		   	fprintf(fCzI, "%e\t", (dz*i+dz/2.0));
   			fprintf(fCzV, "%e\t", (dz*i+dz/2.0));
   			fprintf(fp, "%e\t", (dz*i+dz/2.0));
     			for (j=0; j<Nspecies; j++){
        			fprintf(fCzI,"%e\t", CzI[i + j*Ndistz]);
        			fprintf(fCzV,"%e\t", CzV[i + j*Ndistz]);
        			fprintf(fp,"%e\t", Fnep[i + j*Ndistz]);
      		}
      		fprintf(fCzI,"\n");
      		fprintf(fCzV,"\n");
      		fprintf(fp,"\n");
		 	};
         break;

      case 2:  // Close files
			fclose(fCzion);
			fclose(fQTz);
			fclose(fCzI);
			fclose(fCzV);
			fclose(fp);
         break;
   }  // end of: switch (action)

} // end of: saveCzdata

void  writelogfile(double Rpmean, double  Rpathmean,
                   double totbackscatt, double tottransmitt,
                   double totsputter, double Ebacktrans,
                   unsigned int totcol, unsigned int iontotcol,
                   unsigned int simcollev, unsigned int totsimcol,
                   unsigned int nofound,unsigned int  totNions)
{
	FILE *fid;
   int i;

   if ((fid = fopen("siimplLog.txt", "w")) == NULL) {
   	printf("Can not open: simLog.txt for writing\n");
      exit(1);
   };
   fprintf(fid,"SIIMPL, simulation of ion implantation, (c) M.S. Janson (2003)\n");
   fprintf(fid,"Log file:\n\n");

	fprintf(fid,"Z2");
   for (i=1; i <= Nspecies; i++)
   	fprintf(fid,"  %d",atom[i].Z);
   fprintf(fid,"\n");

	fprintf(fid,"M2");
   for (i=1; i <= Nspecies; i++)
   	fprintf(fid,"  %f",atom[i].M);
   fprintf(fid,"  (amu)\n\n");

 	fprintf(fid, "Number of sub cells (SC) = %d\n", NSC);
   fprintf(fid, "Totat number of atoms in all SC = %d\n", NatSC);
   fprintf(fid, "Avarege number of atoms per SC = %f\n\n", (double)NatSC/NSC);

   fprintf(fid,"Z1 %d\nM1 %f /(amu)\nE0 %f (keV)\ndose %e(cm^-2)\n\n",
   			atom[0].Z, atom[0].M, ion0.E, dose);

	fprintf(fid,"Simulated pseudo ions        : %d\n",Nions);
 	fprintf(fid,"  (including rare event ions): %d\n\n",totNions);

   fprintf(fid,"%% Primary ion ESP data:\n");
   fprintf(fid,"esp.A1 %f (eV/(10^15cm^-2))\n",esp[0].A1*1e2);
   fprintf(fid,"esp.A2 %f \t(= esp.p)\n",esp[0].p);
   fprintf(fid,"esp.p2 %f\n", esp[0].p2);
   fprintf(fid,"esp.A3 %f (keVeV/(10^15cm^-2))\n",esp[0].A3*1e-1);
   fprintf(fid,"esp.A4 %f (keV^-1)\n",esp[0].A4*1e-3);
   fprintf(fid,"esp.A5 %f (keV)\n\n",esp[0].A5*1e3);

	fprintf(fid,"esp.fl %f\nesp.s  %f\n\n",
   				 esp[0].lf, esp[0].s);


   fprintf(fid,"IAPfunction: %d\n",IAPfunction);


   fprintf(fid,"pmax = %f (Å)\n",pmax);
   fprintf(fid,"Estop = %e (keV)\n",Estop);
   fprintf(fid,"vibalg: %d\n\n",vibalg);

   fprintf(fid,"Implanted direction:\n");
   fprintf(fid,"%e\t%e\t%e\n\n",ion0.ev0.x,ion0.ev0.y,ion0.ev0.z);

   fprintf(fid,"Surface Normal:\n");
   fprintf(fid,"%e\t%e\t%e\n\n",ensurf.x,ensurf.y,ensurf.z);

	fprintf(fid,"Implanted rectangle area (Å):\n");
   fprintf(fid,"%e\t%e\t%e\n",Rimp.x,Rimp.y,Rimp.z);
   fprintf(fid,"%e\t%e\t%e\n", Rimp.x + impArea1*e1surf.x,
   	Rimp.y + impArea1*e1surf.y,	Rimp.z + impArea1*e1surf.z);
   fprintf(fid,"%e\t%e\t%e\n",
   	Rimp.x + impArea1*e1surf.x + impArea2*e2surf.x,
   	Rimp.y + impArea1*e1surf.y + impArea2*e2surf.y,
   	Rimp.z + impArea1*e1surf.z + impArea2*e2surf.z);
   fprintf(fid,"%e\t%e\t%e\n\n", Rimp.x + impArea2*e2surf.x,
   	Rimp.y + impArea2*e2surf.y,	Rimp.z + impArea2*e2surf.z);

   fprintf(fid,"surfaceEnergy = %f (keV)\n\n", U_s);


   if (simcollalg) {
   	fprintf(fid,"Simultaneous collision algorithm with dRsc = %f\n",dRsc);
	   fprintf(fid,"Total collisions in simultaneus collisions: %d\n",totsimcol);
      if (simcollev != 0 )
   		fprintf(fid,"Avarege collisions in simultaneous collision event: %e\n",
   			(float) totsimcol/simcollev);
   }

	if (recoil)
		fprintf(fid,"Damage model: Full cascade\n\n");
	else
		fprintf(fid,"Damage model: Kinchin-Pease\n\n");

	if (!randomtarget)
	{
	   fprintf(fid,"Random surface layer:%f Å\n\n",randsurflayer);
	 	fprintf(fid,"rndscatfac: %e\n",rndscatfac);
   }
   else
   	fprintf(fid,"Random ordered, amorphous target\n\n");

   fprintf(fid,"Total number of collisions: %d\n",totcol);
   fprintf(fid,"Total number of collisions for primary ion: "
   	"%d\n",iontotcol);

   fprintf(fid,"total number of nofound: %d\n\n",nofound);

   fprintf(fid,"Rp %f\nRpath %f\n\n", Rpmean, Rpathmean);

	if (NrareDepth == 0)
      fprintf(fid,"<dRv0> (ion) =  %f\n", (Rpathmean*Nions/iontotcol) );
   else
   	fprintf(fid,"<dRv0> not calculated with rare depth algorithm.\n" );


   fprintf(fid,"Random dRv0 =  %f\n\n", (1.0/(nAtom*1e-24)/(pi*sqr(pmax))) );

   fprintf(fid,"Backscattered ions: %e (cm^-2)\n",totbackscatt);
   fprintf(fid,"Transmitted ions: %e (cm^-2)\n",tottransmitt);
   fprintf(fid,"Sputtered target atoms: %e (cm^-2)\n\n",totsputter);

   fprintf(fid,"Energy test: Eout/Ein - 1 = %e\n\n",
                 (energytest(Ebacktrans)/ion0.E) - 1.0 );



   fclose(fid);
}


//  R3 files functions
void	openionR3file(void)
{
	if ((fionR3 = fopen("ionR3.txt", "w")) == NULL) {
   	printf("Can not open: %s for writing\n","ionR3.txt");
      exit(1);
   }
   else
   	fprintf(fionR3, "%% R0.x\tR0.y\tR0.z\tRstop.x\tRstop.y\tRstop.z\n");
}
void	openIVR3files(void)
{
   	if ((fIR3=fopen("IR3.txt", "w")) == NULL) {
   		printf("Can not open: %s for writing\n","IR3.txt");
	      exit(1);
      }
      else
      	fprintf(fIR3, "%\tI.id\tRstop.x\tRstop.y\tRstop.z\n");

   	if ((fVR3=fopen("VR3.txt", "w")) == NULL) {
   		printf("Can not open: %s for writing\n","VR3.txt");
	      exit(1);
      }
      else
      	fprintf(fVR3, "%\tV.id\tR.x\tR.y\tR.z\n");
}
// -------------------------------------------------------




// Help functions  -------------------------------------
void fskipline(FILE *fid)
{
	int ch;

   do {
      ch = getc(fid);
   } while ((!feof(fid)) && ( !((ch==10) || (ch==13)) ))  ;
};
//------------------------------------------------------


