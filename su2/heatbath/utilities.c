#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "algebra.h"
#endif

int N;

double * getelv(double * vector , int i){
	return vector + i;
}

double complex getelm(double * vector , int i, int j){
	if(i==0){
		if(j==0){
			return *getelv(vector , 0) + I*(*getelv(vector , 3));
		}
		else{
			return *getelv(vector , 2) + I*(*getelv(vector , 1));
		}
	}
	else{
		if(j==0){
			return -*getelv(vector , 2) + I*(*getelv(vector , 1));
		}
		else{
			return *getelv(vector , 0) - I*(*getelv(vector , 3));
		}
	}
}

void printv(double * vector){
	printf("(%lf , %lf , %lf , %lf) \n", vector[0],vector[1],vector[2],vector[3]);
}

void printr(double r){
	printf(" %e ", r ); 
}		

void printc(double complex z){
	//function to print a complex doubles
	printf(" %lf + %lf i ", creal(z) , cimag(z) ); 
}			  

void fprintc( FILE * file , double complex z){
	fprintf(file, " %lf + %lf i", creal(z) , cimag(z) );
}

void printm(double * A){
	//function to print a square matrix (2x2) of complex doubles
	int i,j;
	for(i=0 ; i<2 ; i++){
		for(j=0 ; j<2 ; j++){
			printc( getelm(A,i,j) );
		}
		printf("\n");
	}
	printf("\n");
}

void fprintm(FILE * file , double * A){
	int i,j;
	for(i=0 ; i<2 ; i++){
		for(j=0 ; j<2 ; j++){
			fprintc(file,  getelm(A,i,j) );
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

double * getU(double * lattice, int t , int x , int y, int z , int mi){
	return ( lattice + t*N*N*N*4*4 + x*N*N*4*4 + y*N*4*4 + z*4*4 + mi*4);
}

void printLattice(double * lattice){
	int t,x,y,z,mi;
	printf("------------------");
	
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						printf("\n");
						printm( getU(lattice,t,x,y,z,mi) );
					}
				}
			}
		}
	}
	printf("------------------");
}

void copyv(double * out , double * in){
	//copy the values of the elements of in matrix to the out matrix
	int i;
	for(i=0 ; i<4 ; i++){
		*( getelv(out,i) ) = *( getelv(in,i) );
	}
}

double * randSU2v(double * vector ,long * seed , double parameter){
	int j;
	float eps = parameter-0.01;				//displacement from identity
	double r[3] , x[4];	
		
	x[0] = signal( ran0(seed) - 0.5 )*sqrt(1-eps*eps);

	for(j=0 ; j < 3 ; j++){
		r[j] = ran0(seed) - 0.5;
	}

	for(j=1 ; j < 4 ; j++){
		x[j] = eps*r[j-1]/sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
	}

	for(j=0 ; j < 4 ; j++){
		*getelv(vector , j) = x[j];
	}	
}

void initl(double * lattice , float parameter , long * seed){
	int t,x,y,z,mi;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						if(parameter == 0){				
							setidv( getU(lattice, t,x,y,z,mi) );							
						}
						else{							
							randSU2v(getU(lattice, t,x,y,z,mi) , seed , parameter);							
						}
					}
				}
			}
		}
	}
}

void reunitl(double * lattice){
	int t,x,y,z,mi,i;	
	double * aux_link;
	double norm;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						aux_link = getU(lattice, t,x,y,z,mi);						
						norm = sqrt(detv(aux_link));
						//printc(detv(aux_link));
						for(i=0;i<4;i++){			
							*getelv(aux_link , i) /= norm;
						}
						//printr(detv(aux_link));
						//printf("\n");
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

void copyl(double * Lout , double * Lin){
	int t,x,y,z,mi,a,b;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						copyv( getU(Lout, t,x,y,z,mi) , getU(Lin,t,x,y,z,mi) );
						//for(a=0 ; a<4 ; a++){
						//	*(Lout + t*N*N*N*4*4 + x*N*N*4*4 + y*N*4*4 + z*4*4 + mi*4 + a)
						//			= *(Lin + t*N*N*N*4*4 + x*N*N*4*4 + y*N*4*4 + z*4*4 + mi*4 + a);							
						//}
					}
				}
			}
		}
	}
}

double * getg(double * g, int t , int x , int y, int z){
	return ( g + t*N*N*N*4 + x*N*N*4 + y*N*4 + z*4);
}

void printg(double * g){
	int t,x,y,z;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					printm(getg(g,t,x,y,z));
				}
			}
		}
	}
}

long * tseed(long * seed){
	//returns a seed based in the system time		
	*seed = time(0);	
}



