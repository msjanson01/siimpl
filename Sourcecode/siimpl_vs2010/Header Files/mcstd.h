// file: mcstd.h
// Martin Janson, 2000-03-01



// function prototypes;
void r2ucsc(struct ionStruct *ion, struct int3D *UCi);

void updatetgatoml_mjvib(struct double3D *tgatoml, int *Ntgatom, int *SCp1,
             const struct ionStruct *ion);


void	findbc_mjvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             const struct double3D *tgatomlist, const int *Ntgatom,
             const struct ionStruct *ion);

void 	findbc_sc_mjvib(struct double3D *R2UC, struct double3D  *RcUC,
	  		    double *dRv0, double *p, int *tgatomi,
             struct scStruct *sclist, int *Nsimcoll,
             const struct ionStruct *ion);

void	findrandomtg(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             const struct ionStruct *ion, double *n);

void	findbc_robvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             const struct ionStruct *ion);

void	findbc_sc_robvib(struct double3D *R2UC, struct double3D *RcUC,
				 double *dRv0, double *p, int *tgatomi,
             struct scStruct *sclist, int *Nsimcoll,
             const struct ionStruct *ion);


void rareEvent(struct ionStruct *ion);

void mjvibrndnumbers(struct vibStruct *vib);

void Xnorm3D(struct double3D *x3D);

void bcScatter(struct double3D *ev1, double *T, struct double3D *ev2,
               const  struct ionStruct *ion, const struct double3D *R2,
               double M2, const struct double3D *Rcol, double THETA);


double energytest(double Erest);

double randmj(void);

void newrandomseed(void);


