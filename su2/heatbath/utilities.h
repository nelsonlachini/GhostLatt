void printv(double * );

void printr(double );	

void printc(double complex );		  

void fprintc( FILE *  , double complex );

void printm(double * );

void fprintm(FILE *  , double * );

double * getU(double * , int , int , int , int  , int );

void printLattice(double * );

void copyv(double *  , double * );

double * randSU2v(double * vector ,long * seed , double parameter);

void initl(double * , float , long *);

void reunitl(double *);

int getStep( int , int );

void copyl(double *  , double * );

/////////////////////////////////////////////////

long * genSeed();
