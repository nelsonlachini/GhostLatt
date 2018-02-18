//Conventions:
//.SU2 vector refers to the real vector u = (u0,u1,u2,u3);
//.SU2 matrix refers to the matrix A composed from the SU2 vector: A = u0.1 + u1.sx + u2.sy + u3.sz,
//	where sx,sy,sz are the 2x2 sigma matrix and 1 is the 2x2 identity matrix.

double tracev( double * pvector );
	//Returns the trace of a SU2 vector(always real)

void sumv( double * pout , double * B , double * C );
	//Computes the SU2 vector sum B+C and stores it in pout

void mmprodv( double * pout , double  * B, double  * C );
	//Computes the matrix product of two SU2 vectors B.C and stores it in pout


void cmprodv( double * pout , double complex k, double * C );
	//Computes the product of a scalar by the SU2 vector k.C and stores it in pout


void hermcv(double * Adagger , double * A);
	//Computes the hermitian conjugate of the SU2 vector A and stores it in Adagger


double detv(double * A);
	//Returns the determinant of the SU2 vector A (always real)


void setidv(double * A);
	//Sets the SU2 vector A equals to the identity matrix

void setzerov(double * A);
	//Set the SU2 vector A equals to the zero matrix


int signal(double number);
	//Returns -1 if number is negative and +1 else