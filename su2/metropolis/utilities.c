//code to various useful methods
#include <complex.h>
#include <time.h>
#include "algebra.c"

int N;
double complex * lattice;


int signal(double number){
	//returns +1 if element is (zero) positive and -1 else
	if(number < 0)	return -1;	
	return 1;
}

void copyMatrixTo(double complex * out , double complex * in){
	//copy the values of the elements of in matrix to the out matrix
	int i,j;
	for(i=0 ; i<D ; i++){
		for(j=0 ; j<D ; j++){
			*( getElementFromSqMatrix(out,i,j) ) = *( getElementFromSqMatrix(in,i,j) );
		}
	}
}

double complex * getFromLattice(int t , int x , int y, int z , int mi){
	return ( lattice + t*N*N*N*SP_DIM*D*D + x*N*N*SP_DIM*D*D + y*N*SP_DIM*D*D + z*SP_DIM*D*D + mi*D*D);
}

double complex * genRandomSU2(long * seed){
	int j,k;
	float eps = 0.1;
	double x0, r[3] , x[3];
	double complex * matrix;
	
	matrix = malloc( sizeof(double complex)*D*D );
	x0 = signal( ran0(seed) - 0.5 )*sqrt(1-eps*eps);
			
	for(j=0 ; j < 3 ; j++){	
		r[j] = ran0(seed) - 0.5;		
	}

	for(j=0 ; j < 3 ; j++){
		x[j] = eps*r[j]/sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
	}

	for(j=0 ; j < 2 ; j++){
		for(k=0 ; k < 2 ; k++){
			*getElementFromSqMatrix(matrix , j , k) = x0*id2[j][k] + I*x[0]*sx[j][k] + I*x[1]*sy[j][k] + I*x[2]*sz[j][k];				
		}
	}
	
	return matrix;
}

void initLattice(float parameter , long * seed){
	int t,x,y,z,mi,a,b;
	double complex * aux_link;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<SP_DIM; mi++){
						if(parameter == 0){				
							copyMatrixTo( getFromLattice(t,x,y,z,mi) , pid2);
						}
						else{
							aux_link = genRandomSU2(seed);
							copyMatrixTo( getFromLattice(t,x,y,z,mi) , aux_link);
							free(aux_link);		//attention to the gen* functions!!!
						}
					}
				}
			}
		}
	}
}

void reunitLattice(int i){
	int t,x,y,z,mi;
	double complex a1,a2;
	double complex norm;
	double complex * aux_link;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<SP_DIM; mi++){
						aux_link = getFromLattice(t,x,y,z,mi);
						a1 = *(aux_link);
						a2 = *(aux_link +1);
						norm = sqrt( a1*conj(a1) + a2*conj(a2) );
						
						a1 /= norm;
						a2 /= norm;
						*(aux_link) = a1;
						*(aux_link +1) = a2;
						*(aux_link +2) = -conj(a2);
						*(aux_link +3) = conj(a1);																		
					}
				}
			}
		}
	}
}

int getStep( int coord , int dcoord ){		
	//return the lattice coordinate given the increment(useful when there is a chance to get a negative coordinate)
	int aux = (coord + dcoord)%N;
	if(aux<0)	
		return N+aux;
	return aux;
}

void printC(double complex z){
	//function to print a complex doubles
	printf(" %lf + %lf i ", creal(z) , cimag(z) ); 
}			  

void fprintC( FILE * file , double complex z){
	fprintf(file, " %lf + %lf i", creal(z) , cimag(z) );
}

void printMatrix(double complex * A ){
	//function to print a square matrix (DxD) of complex doubles
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

long * genSeed(){
	//returns a seed based in the system time
	long * seed;
	seed = malloc(sizeof(long));
	*seed = time(0);	
	return seed;
}
