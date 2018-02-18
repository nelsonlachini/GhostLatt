//Notations:
//.SU2 vector refers to the real vector u = (u0,u1,u2,u3);
//.SU2 matrix refers to the matrix A composed from the SU2 vector: A = u0.1 + u1.sx + u2.sy + u3.sz,
//where sx,sy,sz are the 2x2 sigma matrix and 1 is the 2x2 identity matrix.

double * getelv(double * A , int i);
	//Returns the i element of the SU2 vector A.

double complex getelm(double * vector , int i, int j);
	//Returns the (i,j) element of the SU2 vector A composed as a SU2 matrix.

void printv(double * A);
	//Prints the SU2 vector A as a vector.

void printr(double r);
	//Prints the double r in exponential notation.

void printc(double complex z);
	//Prints the complete double complex z.

void fprintc( FILE * f , double complex z);
	//Writes the complete double complex z in file f.

void printm(double * A);
	//Prints the SU2 vector A as a SU2 matrix.

void fprintm(FILE * f , double * A);
	//Writes the SU2 vector A as a SU2 matrix in file f.

double * getU(double * L, int t , int x , int y, int z , int mu);
	//Returns the SU2 vector U_mu(t,x,y,z) (a link variable) from the lattice L={U}.

void printLattice(double * L);
	//Prints the entire lattice L as SU2 matrices.

void copyv(double * vout , double * vin);
	//Copies the SU2 vector 'vin' into the SU2 vector 'vout'.

void randSU2v(double * A ,long * seed , double eps);
	//Generates a SU2 vector and stores it in A;
	//seed is the seed pointer used in the random number generator;
	//eps takes any real value inside the interval [0,1]
	//	eps = 0 corresponds to the generation of an identity SU2 vector
	//	eps = 1 corresponds to the random generation of an SU2 vector as far as possible from the identity, i.e. u0=0	

void initl(double * L , float eps , long * seed);
	//Initializes the lattice L with the 'temperature' eps (see randSU2v for details);
	//seed is the seed pointer used in the random number generator.

void reunitl(double * L);
	//Reunitarizes the entire lattice L, i.e. enforces that all its links have determinant equal to 1.

int getStep( int coord , int dcoord );
	//Returns the lattice coordinate given the current coordinate coord and the increment dcoord
	//	(useful when there is a chance to get a negative coordinate).

void copyl(double * Lout , double * Lin);
	//Copies the lattice Lin into the lattice Lout.

double * getg(double * G, int t , int x , int y, int z);
	//Returns the SU2 vector g(t,x,y,z) (gauge transformation element) from the gauge transformation G = {g}.

void printG(double * G);
	//Prints the entire gauge transformation G = {g} as SU2 matrices.

long * tseed(long * seed);
	//Returns a seed based in the system time		