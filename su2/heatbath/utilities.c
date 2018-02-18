#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "algebra.h"
#endif

int N;

double * getelv(double * A , int i){
	return A + i;
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

void printv(double * A){
	printf("(%lf , %lf , %lf , %lf) \n", A[0],A[1],A[2],A[3]);
}

void printr(double r){
	printf(" %e ", r ); 
}		

void printc(double complex z){
	printf(" %lf + %lf i ", creal(z) , cimag(z) ); 
}			  

void fprintc( FILE * f , double complex z){
	fprintf(f, " %lf + %lf i", creal(z) , cimag(z) );
}

void printm(double * A){
	int i,j;
	for(i=0 ; i<2 ; i++){
		for(j=0 ; j<2 ; j++){
			printc( getelm(A,i,j) );
		}
		printf("\n");
	}
	printf("\n");
}

void fprintm(FILE * f , double * A){
	int i,j;
	for(i=0 ; i<2 ; i++){
		for(j=0 ; j<2 ; j++){
			fprintc(f,  getelm(A,i,j) );
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
}

double * getU(double * L, int t , int x , int y, int z , int mu){
	return ( L + t*N*N*N*4*4 + x*N*N*4*4 + y*N*4*4 + z*4*4 + mu*4);
}

void printLattice(double * L){
	int t,x,y,z,mi;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						printf("\n");
						printm( getU(L,t,x,y,z,mi) );
					}
				}
			}
		}
	}	
}

void copyv(double * vout , double * vin){
	int i;
	for(i=0 ; i<4 ; i++){
		*( getelv(vout,i) ) = *( getelv(vin,i) );
	}
}

void randSU2v(double * A ,long * seed , double eps){
	int j;
	double r[3] , x[4];
	
	x[0] = signal( ran0(seed) - 0.5 )*sqrt(1-eps*eps);
	for(j=0 ; j < 3 ; j++){
		r[j] = ran0(seed) - 0.5;
	}

	for(j=1 ; j < 4 ; j++){
		x[j] = eps*r[j-1]/sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
	}
	for(j=0 ; j < 4 ; j++){
		*getelv(A , j) = x[j];
	}	
}

void initl(double * L , float eps , long * seed){
	int t,x,y,z,mi;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						if(eps == 0){
							setidv( getU(L, t,x,y,z,mi) );
						}
						else{
							randSU2v(getU(L, t,x,y,z,mi) , seed , eps);
						}
					}
				}
			}
		}
	}
}

void reunitl(double * L){
	int t,x,y,z,mi,i;	
	double * aux_link;
	double norm;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						aux_link = getU(L, t,x,y,z,mi);						
						norm = sqrt(detv(aux_link));						
						for(i=0;i<4;i++){			
							*getelv(aux_link , i) /= norm;
						}
					}
				}
			}
		}
	}
}

int getStep( int coord , int dcoord ){		
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
					}
				}
			}
		}
	}
}

double * getg(double * G, int t , int x , int y, int z){
	return ( G + t*N*N*N*4 + x*N*N*4 + y*N*4 + z*4);
}

void printG(double * G){
	int t,x,y,z;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					printm(getg(G,t,x,y,z));
				}
			}
		}
	}
}

long * tseed(long * seed){
	*seed = time(0);	
}



