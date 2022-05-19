// file: siimpl.h
// Martin Janson, 2000-03-13

// function prototypes

void pseudoIonSim(double *Rp,
						int *backscatt, int *transmitt,
					   struct ionStruct *ion, double *Erest,
                  unsigned int *collisions,
                  unsigned int *simcollev, unsigned int *totsimcol,
                  unsigned int *nofound);

void add2ill(struct ionStruct *ion);




