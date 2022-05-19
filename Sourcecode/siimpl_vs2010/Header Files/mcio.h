// file: mcio.h
// Martin Janson, 2000-03-12


// function prototypes;
void fskipline(FILE *fid);

void structureSize(char *filename, int *NatUC,
											  int *Nspecies, int *Ndz);

int readSimData(char *filename);

int readnormdist(char *filename);

void exitRD(char *mesg);
void exitRD2many(char *mesg);

//  ---- Import Cz data functions
void 	importCzion(const char path[100]);

void 	importCzI(const char path[100]);

void 	importCzV(const char path[100]);



// output functions
void saveCzData(int action);

void  writelogfile(double Rpmean, double  Rpathmean,
                   double totbackscatt, double tottransmitt,
                   double totsputter, double Ebacktrans,
                   unsigned int totcol, unsigned int iontotcol,
                   unsigned int simcollev, unsigned int totsimcol,
                   unsigned int nofound, unsigned int  totNions);

void	openionR3file(void);

void	openIVR3files(void);


