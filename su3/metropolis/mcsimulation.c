#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran0.c"
#include "algebra.c"

#define POOL_SIZE 2

int N,N_cf,N_hit;
double avg_acceptance;
double beta,eps;

double complex * lattice;
double complex * pool;

//return the adress of the link in the random pool given the index
double complex * getFromPool(int index){
	return(pool + index*D*D);
}

void printLattice(double complex * mylattice){
	int t,x,y,z,mi;
	double complex * aux_link;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<SP_DIM; mi++){
						aux_link = mylattice + t*N*N*N*SP_DIM*D*D + x*N*N*SP_DIM*D*D + y*N*SP_DIM*D*D + z*SP_DIM*D*D + mi*D*D;
						printf("\n Link ::: %d%d%d%d%d \n",t,x,y,z,mi);
						printMatrix(aux_link);
					}
				}
			}
		}
	}
}

//algorithm from Lang's book for generating random su2 matrices close to the identity, controlled by the parameter eps
void createPoolSU3(long * seed){
	double x0;
	double s[3] , r[3] , t[3] , x[3];
	double complex * R;
	double complex * S;
	double complex * T;
	R = genId();
	S = genId();
	T = genId();
	int i,j,k;
		
	double complex * aux_hermitian_conj;	
	
	for(i=0 ; i < POOL_SIZE ; i++){
		//generating R matrix
		x0 = signal( ran0(seed) - 0.5 )*sqrt(1 - eps*eps);			
		for(j=0 ; j < 3 ; j++){	
			x[j] = ran0(seed) - 0.5;		
		}		
		for(j=0 ; j < 3 ; j++){
			r[j] = eps*x[j]/sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
		}
		*getElementFromSqMatrix(R , 0 , 0) = x0 + I*r[2]*sz[0][0];
		*getElementFromSqMatrix(R , 1 , 1) = x0 + I*r[2]*sz[1][1];
		*getElementFromSqMatrix(R , 0 , 1) = I*r[0]*sx[0][1] + I*r[1]*sy[0][1];
		*getElementFromSqMatrix(R , 1 , 0) = I*r[0]*sx[1][0] + I*r[1]*sy[1][0];
		
		//generating S matrix
		x0 = signal( ran0(seed) - 0.5 )*sqrt(1 - eps*eps);			
		for(j=0 ; j < 3 ; j++){	
			x[j] = ran0(seed) - 0.5;		
		}		
		for(j=0 ; j < 3 ; j++){
			s[j] = eps*x[j]/sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
		}
		*getElementFromSqMatrix(S , 0 , 0) = x0 + I*s[2]*sz[0][0];
		*getElementFromSqMatrix(S , 2 , 2) = x0 + I*s[2]*sz[1][1];
		*getElementFromSqMatrix(S , 0 , 2) = I*s[0]*sx[0][1] + I*s[1]*sy[0][1];
		*getElementFromSqMatrix(S , 2 , 0) = I*s[0]*sx[1][0] + I*s[1]*sy[1][0];
		
		
		//generating T matrix
		x0 = signal( ran0(seed) - 0.5 )*sqrt(1 - eps*eps);			
		for(j=0 ; j < 3 ; j++){	
			x[j] = ran0(seed) - 0.5;		
		}		
		for(j=0 ; j < 3 ; j++){
			t[j] = eps*x[j]/sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
		}
		*getElementFromSqMatrix(T , 1 , 1) = x0 + I*t[2]*sz[0][0];
		*getElementFromSqMatrix(T , 2 , 2) = x0 + I*t[2]*sz[1][1];
		*getElementFromSqMatrix(T , 1 , 2) = I*t[0]*sx[0][1] + I*t[1]*sy[0][1];
		*getElementFromSqMatrix(T , 2 , 1) = I*t[0]*sx[1][0] + I*t[1]*sy[1][0];
				
		matrixProductInout(getFromPool(i) , R , S);
		matrixProductInout(getFromPool(i) , getFromPool(i) , T);		
		aux_hermitian_conj = genHermitianConj( getFromPool(i) );
		copyMatrixTo( getFromPool(i+POOL_SIZE) , aux_hermitian_conj );
		//free(aux_hermitian_conj);		//always free the memory after using gen* methods!!!
	}
	free(R);
	free(S);
	free(T);
	
}

void createPoolElton(long * seed){
	double U0, U1, U2, U3, ro , phi;
	double U[3], x[3];
	int i,j,k;
	
	double complex * aux_hermitian_conj;
	
	for(i=0 ; i < POOL_SIZE ; i++){
		U0 = ran0(seed);		
		ro = ran0(seed);
		phi = ran0(seed);
		
		U1 = sqrt(1- ro*ro)*sqrt(1-U0*U0)*cos(phi);
		U2 = sqrt(1- ro*ro)*sqrt(1-U0*U0)*sin(phi);
		U3 = ro*sqrt(1-U0*U0);
		
		for(j=0 ; j < 2 ; j++){
			for(k=0 ; k < 2 ; k++){
				*getElementFromSqMatrix(getFromPool(i) , j , k) = U0*id2[j][k] + I*U1*sx[j][k] + I*U2*sy[j][k] + I*U3*sz[j][k];				
			}
		}
		aux_hermitian_conj = genHermitianConj( getFromPool(i) );
		copyMatrixTo( getFromPool(i+POOL_SIZE) , aux_hermitian_conj );
		free(aux_hermitian_conj);		//always free the memory after using gen* methods!!!
	}	
}

double complex * getFromLattice(int t , int x , int y, int z , int mi){
	return ( lattice + t*N*N*N*SP_DIM*D*D + x*N*N*SP_DIM*D*D + y*N*SP_DIM*D*D + z*SP_DIM*D*D + mi*D*D);
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

void initLattice(float parameter){
	int t,x,y,z,mi,a,b;
	double complex * aux_link;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<SP_DIM; mi++){
						aux_link = getFromLattice(t,x,y,z,mi);
						for(a=0 ; a<D ; a++){				//cold start
							for(b=0 ; b<D ; b++){								
								*(aux_link + D*a + b) = (float)(a==b);
							}
						}						
					}
				}
			}
		}
	}
}

double complex calcTrace( double complex * pmatrix ){
	int i;
	double complex out = 0;
	for(i=0 ; i < D ; i++){		
		out += *( pmatrix + i*D + i ); // (i,i) element of pmatrix
	}
	return out;
}

double calcAction( double complex * link , double complex * _staple ){
	double ans;
	double complex * product = genMatrixProduct( link , _staple );	
	ans = creal( calcTrace( product ));
	free(product);				//always use free after using gen* methods!!!!!!!!!!!!!!!!
	
	return -beta*ans/(float)(D);
}

int getStep( int coord , int dcoord ){		//return the lattice coordinate given the increment(useful when there is a chance to get a negative coordinate)
		int aux = (coord + dcoord)%N;
		if(aux<0)	
			return N+aux;
		return aux;
}

double updateStaple(int t, int x, int y, int z, int mi , double complex * _staple){
	int i,ni;
	int ami[SP_DIM] = {0,0,0,0};
	int ani[SP_DIM];
    ami[mi] = 1;
	
	double complex * auxp = malloc((sizeof(double complex))*D*D);
	double complex * aux1 ;
	double complex * aux2 ;
	double complex * aux3;
	
	copyMatrixTo( _staple , pzero3 ); 		//now staple should be a zero DxD matrix
	
    for(ni=0 ; ni<SP_DIM ; ni++){
        if(ni != mi){
			for(i=0 ; i<SP_DIM ; i++)	ani[i] = 0;	
            ani[ni] = 1 ;
			
			aux1 = getFromLattice( (t+ami[0])%N , (x+ami[1])%N , (y+ami[2])%N , (z+ami[3])%N, ni) ;
			aux2 = genHermitianConj( getFromLattice( (t+ani[0])%N , (x+ani[1])%N , (y+ani[2])%N , (z+ani[3])%N , mi));
			aux3 = genHermitianConj( getFromLattice(t , x , y , z , ni) );
	
			matrixProductInout(auxp, aux1 , aux2 );
			matrixProductInout(auxp, auxp , aux3);
			
			free(aux2);
			free(aux3);
			
			matrixSumInout(_staple , _staple , auxp);
            					
			aux1 = genHermitianConj(getFromLattice( getStep(t,ami[0]-ani[0]), getStep(x,ami[1]-ani[1]) , getStep(y,ami[2]-ani[2]) , getStep(z,ami[3]-ani[3]) , ni));
			aux2 = genHermitianConj(getFromLattice( getStep(t,-ani[0]) , getStep(x,-ani[1]) , getStep(y,-ani[2]) , getStep(z,-ani[3]) , mi));
			aux3 = getFromLattice( getStep(t,-ani[0]) , getStep(x,-ani[1]) , getStep(y,-ani[2]) , getStep(z,-ani[3]) , ni);
														
			matrixProductInout(auxp , aux1 , aux2);
			matrixProductInout(auxp , auxp , aux3);
			
			free(aux1);
			free(aux2);
			
			matrixSumInout(_staple , _staple , auxp);			
		}
	}
	free(auxp);;
}

double complex * getRandomFromPool(long * _seed){		
	return getFromPool((int)roundf(ran0(_seed)*(2*POOL_SIZE-1)));
}

void updateLink(int t, int x, int y, int z, int mi , long * seed){
	int i;
	double old_S , new_S , dS;
	double complex * potential_link;
	double complex * old_link = getFromLattice(t,x,y,z,mi);
	double complex * staple = malloc( sizeof(double complex)*D*D );
	double complex * aux_link;
	updateStaple(t,x,y,z,mi , staple);
	
	for(i=0 ; i < N_hit ; i++){
		old_S = calcAction( old_link , staple );
		potential_link = genMatrixProduct( getRandomFromPool(seed) , old_link);
		new_S = calcAction( potential_link , staple);
		
		dS = new_S - old_S;
		if(  ran0(seed) < exp(-dS)  ){        //case where the update is approved
			copyMatrixTo( old_link , potential_link);
		}
		free(potential_link);		//free pointers after gen* calls
	}
	free(staple);
}

void updateLattice(long * seed){	
	int t,x,y,z,mi;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<SP_DIM ; mi++){
						updateLink(t,x,y,z,mi , seed );						
					}
				}
			}
		}
	}
}

double wilsonLoopMeasure(int X, int Y ){
    double sum = 0;
	int temp;
	int k,j,mi,ni,t,x,y,z;
    int ani[4],ami[4];
	double complex * auxc;
	double complex * aux = malloc(sizeof(double complex)*D*D);
    for(t=0 ; t<N ; t++){
        for(x=0 ; x<N ; x++){    
            for(y=0 ; y<N ;y++){
                for(z=0 ; z<N ; z++){                    
                    for(ni=0 ; ni<3 ; ni++){
                        for(k=0 ; k<SP_DIM ; k++)	
							ani[k] = 0;
                        ani[ni] = 1;
                        for(mi=ni+1 ; mi<4 ; mi++){ 
                            for(k=0 ; k<4 ; k++)	
								ami[k] = 0;
                            ami[mi] = 1;                            
                            k=0;
                            j=0;                         
                            
                            copyMatrixTo( aux , pid2 ); 
                                                                                 
                            for(k=0 ; k<X ; k++){
								matrixProductInout(aux , aux ,
									getFromLattice( (t+k*ami[0])%N , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
							}
	
							k--;
                            for(j=0 ; j<Y ; j++){
								matrixProductInout(aux , aux , getFromLattice( (t+(k+1)*ami[0] + j*ani[0])%N , (x+(k+1)*ami[1]+ j*ani[1])%N ,
										(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
							}
	
							j--;
							for( ; k > -1 ; k--){
								auxc = genHermitianConj( getFromLattice( (t+k*ami[0] + (j+1)*ani[0])%N , (x+k*ami[1]+ (j+1)*ani[1])%N , 
									(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
								matrixProductInout(aux, aux , auxc);
								free(auxc);
							}
	
							k++;
							for( ; j > -1 ; j--){
								auxc = genHermitianConj(getFromLattice( (t+ j*ani[0])%N , (x+ j*ani[1])%N , 
									(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
								matrixProductInout(aux , aux , auxc);
								free(auxc);
							}

							sum += creal(calcTrace(aux))/D;							
						}
					}
				}
			}
		}
	}
	
	if (X != Y){
		temp = X;
		X = Y;
		Y = temp;
		for(t=0 ; t<N ; t++){
			for(x=0 ; x<N ; x++){    
				for(y=0 ; y<N ;y++){
					for(z=0 ; z<N ; z++){                    
						for(ni=0 ; ni<3 ; ni++){
							for(k=0 ; k<SP_DIM ; k++)	
								ani[k] = 0;
							ani[ni] = 1;
							for(mi=ni+1 ; mi<4 ; mi++){ 
								for(k=0 ; k<4 ; k++)	
									ami[k] = 0;
								ami[mi] = 1;                            
								k=0;
								j=0;                         

								copyMatrixTo( aux , pid2 ); 

								for(k=0 ; k<X ; k++){
									matrixProductInout(aux , aux ,
										getFromLattice( (t+k*ami[0])%N , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
								}

								k--;
								for(j=0 ; j<Y ; j++){
									matrixProductInout(aux , aux , getFromLattice( (t+(k+1)*ami[0] + j*ani[0])%N , (x+(k+1)*ami[1]+ j*ani[1])%N ,
											(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
								}

								j--;
								for( ; k > -1 ; k--){
									auxc = genHermitianConj( getFromLattice( (t+k*ami[0] + (j+1)*ani[0])%N , (x+k*ami[1]+ (j+1)*ani[1])%N , 
										(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
									matrixProductInout(aux, aux , auxc);
									free(auxc);
								}

								k++;
								for( ; j > -1 ; j--){
									auxc = genHermitianConj(getFromLattice( (t+ j*ani[0])%N , (x+ j*ani[1])%N , 
										(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
									matrixProductInout(aux , aux , auxc);
									free(auxc);
								}

								sum += creal(calcTrace(aux))/D;							
							}
						}
					}
				}
			}
		}
		free(aux);
		return sum/(12*N*N*N*N);
	}
    
	free(aux);
	return sum/(6*N*N*N*N);
}