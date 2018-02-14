//code to various useful methods
#include <complex.h>
#include <time.h>

#define D 3
#define SP_DIM 4
int N;

//function to print a complex doubles
void printC(double complex z){
	printf(" %lf + %lf i ", creal(z) , cimag(z) ); 
}			  

void fprintC( FILE * file , double complex z){
	fprintf(file, " %lf + %lf i", creal(z) , cimag(z) );
}


//return the adress of the D dimensional square matrix element indicated
double complex * getElementFromSqMatrix(double complex * matrix , int i , int j){
	return(matrix + i*D + j);
}

double complex * getElementFromSqMatrix2(double complex * matrix , int i , int j){
	return(matrix + i*2 + j);
}

//copy the values of the elements of in matrix to the out matrix
void copyMatrixTo(double complex * out , double complex * in){
	int i,j;
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*( getElementFromSqMatrix(out,i,j) ) = *( getElementFromSqMatrix(in,i,j) );
		}
	}
}

void copyLatticeTo(double complex * Lout , double complex * Lin){
	int t,x,y,z,mi,a,b;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<SP_DIM; mi++){
						for(a=0 ; a<D ; a++){
							for(b=0 ; b<D ; b++){
								*(Lout + t*N*N*N*SP_DIM*D*D + x*N*N*SP_DIM*D*D + y*N*SP_DIM*D*D + z*SP_DIM*D*D + mi*D*D + a*D +b)
									= *(Lin + t*N*N*N*SP_DIM*D*D + x*N*N*SP_DIM*D*D + y*N*SP_DIM*D*D + z*SP_DIM*D*D + mi*D*D + a*D +b);
							}
						}
					}
				}
			}
		}
	}
}


//returns a seed based in the system time
long * genSeed(){
	long * seed;
	seed = malloc(sizeof(long));
	*seed = time(0);	
	return seed;
}
