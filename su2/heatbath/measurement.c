#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utilities.h"

#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "algebra.h"
#endif 

int N;
double beta;

double wilsonLoopMeasure(double * lattice , int X, int Y ){
    double sum = 0;
	int temp;
	int k,j,mi,ni,t,x,y,z;
    int ani[4],ami[4];
	double * auxc = malloc(sizeof(double)*4);
	double * aux = malloc(sizeof(double)*4);
    for(t=0 ; t<N ; t++){
        for(x=0 ; x<N ; x++){    
            for(y=0 ; y<N ;y++){
                for(z=0 ; z<N ; z++){                    
                    for(ni=0 ; ni<3 ; ni++){
                        for(k=0 ; k<4 ; k++) ani[k] = 0;
                        ani[ni] = 1;
                        for(mi=ni+1 ; mi<4 ; mi++){ 
                            for(k=0 ; k<4 ; k++) ami[k] = 0;
                            ami[mi] = 1;                            
                            k=0;
                            j=0;                         
                            
                            setidv( aux ); 
                                                                                 
                            for(k=0 ; k<X ; k++){
								mmprodv(aux , aux ,
									getU(lattice, (t+k*ami[0])%N , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
							}
	
							k--;
                            for(j=0 ; j<Y ; j++){
								mmprodv(aux , aux , getU(lattice, (t+(k+1)*ami[0] + j*ani[0])%N , (x+(k+1)*ami[1]+ j*ani[1])%N ,
										(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
							}
	
							j--;
							for( ; k > -1 ; k--){
								hermcv( auxc , getU(lattice, (t+k*ami[0] + (j+1)*ani[0])%N , (x+k*ami[1]+ (j+1)*ani[1])%N , 
									(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
								mmprodv(aux, aux , auxc);								
							}
	
							k++;
							for( ; j > -1 ; j--){
								hermcv( auxc , getU(lattice, (t+ j*ani[0])%N , (x+ j*ani[1])%N , 
									(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
								mmprodv(aux , aux , auxc);								
							}
							sum += tracev(aux)/2.;							
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
							for(k=0 ; k<4 ; k++) ani[k] = 0;
							ani[ni] = 1;
							for(mi=ni+1 ; mi<4 ; mi++){ 
								for(k=0 ; k<4 ; k++) ami[k] = 0;
								ami[mi] = 1;                            
								k=0;
								j=0;                         

								setidv( aux ); 

								for(k=0 ; k<X ; k++){
									mmprodv(aux , aux ,
										getU(lattice, (t+k*ami[0])%N , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
								}

								k--;
								for(j=0 ; j<Y ; j++){
									mmprodv(aux , aux , getU(lattice, (t+(k+1)*ami[0] + j*ani[0])%N , (x+(k+1)*ami[1]+ j*ani[1])%N ,
											(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
								}

								j--;
								for( ; k > -1 ; k--){
									hermcv( auxc , getU(lattice, (t+k*ami[0] + (j+1)*ani[0])%N , (x+k*ami[1]+ (j+1)*ani[1])%N , 
										(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
									mmprodv(aux, aux , auxc);									
								}

								k++;
								for( ; j > -1 ; j--){
									hermcv( auxc , getU(lattice, (t+ j*ani[0])%N , (x+ j*ani[1])%N , 
										(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
									mmprodv(aux , aux , auxc);									
								}
								sum += tracev(aux)/2.;							
							}
						}
					}
				}
			}
		}
		free(aux);
		free(auxc);
		return sum/(12*N*N*N*N);
	}
    free(auxc);
	free(aux);
	return sum/(6*N*N*N*N);
}