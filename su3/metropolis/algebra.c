#include <complex.h>
#include <stdio.h>
#include "utilities.c"


//code to algebra methods
//Note: in my convention the gen* methods always allocate some memory to operate, being needed to free such memory at sometime after

////////////////////definition of useful CONSTANT matrices
//definition of the 2 dimensional identity matrix according with the dimension D defined
double complex id2[2][2] = {{1,0},{0,1}};
double complex id3[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

double complex * pid2 = &id2[0][0];
double complex * pid3 = &id3[0][0];


//definition of the sigma matrices			  
double complex sx[2][2] = {{0,1},{1,0}};
double complex sy[2][2] = {{0,-I},{I,0}};
double complex sz[2][2] = {{1,0},{0,-1}};

double complex * psx = &sx[0][0];
double complex * psy = &sy[0][0];
double complex * psz = &sz[0][0];

double complex zero2[2][2] = {{0,0},{0,0}};
double complex * pzero2 = &zero2[0][0];
double complex zero3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double complex * pzero3 = &zero3[0][0];

/////////////////////////definition of useful algebra methods
//function to print a square matrix (DxD) of complex doubles
void printMatrix(double complex * A ){
	int i,j;
	double complex element;
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			printC( *(getElementFromSqMatrix(A,i,j)) );
		}
		printf("\n");
	}
	printf("\n");
}

void fprintMatrix(FILE * file , double complex * A){
	int i,j;
	double complex element;
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			fprintC(file,  *(getElementFromSqMatrix(A,i,j)) );
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

//function to sum two matrices of complex double elements, B and C; result goes into *pout (memory location that pout points)
//OBS.: *B and *C are not modified here!
void matrixSumInout( double complex * pout , double complex * B , double complex * C ){
	double complex A[D][D];
	int i,j;		
	
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			A[i][j] = *getElementFromSqMatrix(B,i,j) + *getElementFromSqMatrix(C,i,j);
		}
	}

	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*getElementFromSqMatrix(pout,i,j) = A[i][j];
		}
	}
}

//function to multiply two matrices of complex double elements, B and C; result goes into *pout (memory location that pout points)
//OBS.: *B and *C are not modified here!
void matrixProductInout( double complex * pout , double complex * B, double complex * C ){
	double complex A[D][D];
	int i,j,k;		
	
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			A[i][j] = 0;
			for(k=0 ; k<D ; k++){						
				A[i][j] += ( *(getElementFromSqMatrix(B,i,k)) )*( *(getElementFromSqMatrix(C,k,j)) );
			}
		}
	}
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*(getElementFromSqMatrix(pout,i,j)) = A[i][j];
		}
	}
}

//function to multiply two matrices of complex double elements, B and C; result is returned in allocated memory as a pointer
//OBS.: *B and *C are not modified here!
double complex * genMatrixProduct(double complex * B , double complex * C ){
	double complex * A;
	A = malloc(sizeof(double complex)*D*D);
	int i,j,k;		
	
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*(getElementFromSqMatrix(A,i,j)) = 0;
			for(k=0 ; k<D ; k++){						
				*(getElementFromSqMatrix(A,i,j)) += *(getElementFromSqMatrix(B,i,k)) * ( *(getElementFromSqMatrix(C,k,j)) );
			}
		}
	}
	return A;
}

double complex * genHermitianConj(double complex * A ) {
	int i,j;
	double complex * Adagger;
	Adagger = malloc( sizeof(double complex)*D*D );
	
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*(getElementFromSqMatrix(Adagger,i,j)) = conj( *(getElementFromSqMatrix(A,j,i) ) );
		}
	}
	return Adagger;
}

//returns +1 if element is (zero) positive and -1 else
int signal(double number){
	if(number < 0)	return -1;	
	return 1;
}

double complex det2(double complex * A){
	return *getElementFromSqMatrix2(A,0,0)*(*getElementFromSqMatrix2(A,1,1)) - *getElementFromSqMatrix2(A,0,1)*(*getElementFromSqMatrix2(A,1,0));
}

double complex det3(double complex * A){
	return *getElementFromSqMatrix(A,0,0)*(*getElementFromSqMatrix(A,1,1))*(*getElementFromSqMatrix(A,2,2)) 
		+ *getElementFromSqMatrix(A,0,1)*(*getElementFromSqMatrix(A,1,2))*(*getElementFromSqMatrix(A,2,0))
		+ (*getElementFromSqMatrix(A,0,2))*(*getElementFromSqMatrix(A,1,0))*(*getElementFromSqMatrix(A,2,1))
		- (*getElementFromSqMatrix(A,2,0))*(*getElementFromSqMatrix(A,1,1))*(*getElementFromSqMatrix(A,0,2))
		- (*getElementFromSqMatrix(A,2,1))*(*getElementFromSqMatrix(A,1,2))*(*getElementFromSqMatrix(A,0,0))
		- (*getElementFromSqMatrix(A,2,2))*(*getElementFromSqMatrix(A,1,0))*(*getElementFromSqMatrix(A,0,1));
}

double complex * genZero(){
	double complex * zero;
	int i,j;
	zero = malloc(sizeof(double complex));
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*(zero + i*D +j) = 0.0;
		}
	}
	return zero;
}

double complex * genId(){
	double complex * zero;
	int i,j;
	zero = malloc(sizeof(double complex)*D*D);
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*(zero + i*D +j) = (double)(i ==j);
		}
	}
	return zero;
}
