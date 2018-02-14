#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran0.c"
#include "utilities.c"

#define POOL_SIZE 2000

int N,N_cf,N_hit;
double avg_acceptance;
double beta,eps;
double complex * lattice;
double complex * pool;

double complex * getFromPool(int index){
	//return the adress of the link in the random pool given the index
	return(pool + index*D*D);
}

double complex * getRandomFromPool(long * _seed){	
	return getFromPool((int)roundf(ran0(_seed)*(2*POOL_SIZE-1)));
}

void createPool(long * seed){
	//algorithm from Lang's book for generating random su2 matrices close to the identity, controlled by the parameter eps
	double x0;
	double r[3], x[3];
	int i,j,k;
	
	double complex * aux_hermitian_conj;
	
	for(i=0 ; i < POOL_SIZE ; i++){
		x0 = signal( ran0(seed) - 0.5 )*sqrt(1 - eps*eps);
			
		for(j=0 ; j < 3 ; j++){	
			r[j] = ran0(seed) - 0.5;		
		}
		
		for(j=0 ; j < 3 ; j++){
			x[j] = eps*r[j]/sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
		}
		
		for(j=0 ; j < 2 ; j++){
			for(k=0 ; k < 2 ; k++){
				*getElementFromSqMatrix(getFromPool(i) , j , k) = x0*id2[j][k] + I*x[0]*sx[j][k] + I*x[1]*sy[j][k] + I*x[2]*sz[j][k];				
			}
		}
	
		aux_hermitian_conj = genHermitianConj( getFromPool(i) );
		copyMatrixTo( getFromPool(i+POOL_SIZE) , aux_hermitian_conj );
		free(aux_hermitian_conj);		//always free the memory after using gen* methods!!!
	} 				
}

double calcAction( double complex * link , double complex * _staple ){
	double ans;
	double complex * product = genMatrixProduct( link , _staple ) ;
	ans = creal( calcTrace( product ));
	
	free(product);				//always use free after using gen* methods!!!!!!!!!!!!!!!!
	
	return -beta*ans/D;
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
	
	copyMatrixTo( _staple , pzero2 ); 		//now staple should be a zero DxD matrix
	
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
	free(auxp);
}

void updateLink(int t, int x, int y, int z, int mi , long * seed){
	int i;
	double old_S , new_S , dS;
	double complex * potential_link;		
	double complex * old_link = getFromLattice(t,x,y,z,mi);
	double complex * staple = malloc( sizeof(double complex)*D*D );
	double complex * aux_link;
	updateStaple(t,x,y,z,mi , staple);
					
	int acc_r = 0;
	for(i=0 ; i < N_hit ; i++){	
		old_S = calcAction( old_link , staple );			
		potential_link = genMatrixProduct( getRandomFromPool(seed) , old_link);
		
		new_S = calcAction( potential_link , staple);
		
		dS = new_S - old_S;
		
		if(  ran0(seed) < exp(-dS)  ){        //case where the update is approved
			acc_r++;	
			copyMatrixTo( old_link , potential_link);					
		}
		free(potential_link);		//free pointers after gen* calls		
	}	
	free(staple);	
	avg_acceptance += acc_r;
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
    
	free(aux);
	return sum/(6*N*N*N*N);
}
