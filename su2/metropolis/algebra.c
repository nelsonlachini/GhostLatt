#include <complex.h>
#include <stdio.h>

#define D 2
#define SP_DIM 4

//code to algebra methods
//Note: in my convention the gen* methods always allocate some memory to operate, being needed to free such memory at sometime later

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


/////////////////////////definition of useful algebra methods
double complex * getElementFromSqMatrix(double complex * matrix , int i , int j){
	//return the adress of the D dimensional square matrix element indicated
	return(matrix + i*D + j);
}

double complex calcTrace( double complex * pmatrix ){
	int i;
	double complex out = 0;
	for(i=0 ; i < D ; i++){		
		out += *( pmatrix + i*D + i ); // (i,i) element of pmatrix
	}
	return out;
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

void matrixProductByConstInout( double complex * pout , double complex k, double complex * C ){
	double complex A[D][D];
	int i,j;
	
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){			
			A[i][j] = *(getElementFromSqMatrix(C,i,j))*k ;
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

double complex det2(double complex * A){
	return *getElementFromSqMatrix(A,0,0)*(*getElementFromSqMatrix(A,1,1)) - *getElementFromSqMatrix(A,0,1)*(*getElementFromSqMatrix(A,1,0));
}
