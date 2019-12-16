#include "measurement.h"
#include "global.h"

#include "inverters.h"

double wilsonLoopMeasureSym(double * lattice , int X, int Y){
    double sum = 0;
    int temp;
    int k,j,mi,ni,t,x,y,z;
    int ani[4],ami[4];
    double * auxc = malloc(sizeof(double)*4);
    double * aux = malloc(sizeof(double)*4);
    for(t=0 ; t<Nt ; t++){
        for(x=0 ; x<N ; x++){
            for(y=0 ; y<N ;y++){
                for(z=0 ; z<N ; z++){
                    for(ni=0 ; ni<3 ; ni++){
                        setquadv(&ani[0],ni);
                        for(mi=ni+1 ; mi<4 ; mi++){
                            setquadv(&ami[0],mi);
                            k=0;
                            j=0;
                            setidv(aux);

                            for(k=0 ; k<X ; k++){
								mmprodv(aux, aux, getU(lattice, (t+k*ami[0])%Nt , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
							}
							k--;

                            for(j=0 ; j<Y ; j++){
								mmprodv(aux , aux , getU(lattice, (t+(k+1)*ami[0] + j*ani[0])%Nt , (x+(k+1)*ami[1]+ j*ani[1])%N ,
										(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
							}
							j--;

                            for( ; k > -1 ; k--){
								hermcv( auxc , getU(lattice, (t+k*ami[0] + (j+1)*ani[0])%Nt , (x+k*ami[1]+ (j+1)*ani[1])%N ,
									(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
								mmprodv(aux, aux , auxc);
							}
							k++;

							for( ; j > -1 ; j--){
								hermcv( auxc , getU(lattice, (t+ j*ani[0])%Nt , (x+ j*ani[1])%N ,
									(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
								mmprodv(aux , aux , auxc);
							}

                            sum += tracev(aux);
                        }
                    }
                }
            }
        }
    }
    if(X != Y){
		temp = X;
		X = Y;
		Y = temp;
		for(t=0 ; t<Nt ; t++){
			for(x=0 ; x<N ; x++){
				for(y=0 ; y<N ;y++){
					for(z=0 ; z<N ; z++){
						for(ni=0 ; ni<3 ; ni++){
							for(k=0 ; k<4 ; k++) ani[k] = 0;
							ani[ni] = 1;
							for(mi=ni+1 ; mi<4 ; mi++){
                                setquadv(&ami[0],mi);
								k=0;
								j=0;
								setidv( aux );

								for(k=0 ; k<X ; k++){
									mmprodv(aux , aux ,
										getU(lattice, (t+k*ami[0])%Nt , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
								}

								k--;
								for(j=0 ; j<Y ; j++){
									mmprodv(aux , aux , getU(lattice, (t+(k+1)*ami[0] + j*ani[0])%Nt , (x+(k+1)*ami[1]+ j*ani[1])%N ,
											(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
								}

								j--;
								for( ; k > -1 ; k--){
									hermcv( auxc , getU(lattice, (t+k*ami[0] + (j+1)*ani[0])%Nt , (x+k*ami[1]+ (j+1)*ani[1])%N ,
										(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
									mmprodv(aux, aux , auxc);
								}

								k++;
								for( ; j > -1 ; j--){
									hermcv( auxc , getU(lattice, (t+ j*ani[0])%Nt , (x+ j*ani[1])%N ,
										(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
									mmprodv(aux , aux , auxc);
								}

								sum += tracev(aux);
							}
						}
					}
				}
			}
		}
		free(aux);
		free(auxc);
		return sum/(2*12.0*totalV);
    }
    free(auxc);
	free(aux);
	return sum/(2*6.0*totalV);
}

double wilsonLoopMeasureAsym(double * lattice , int R, int T){
    double sum = 0;
    int temp;
    int k,j,mi,ni,t,x,y,z;
    int ani[4],ami[4];
    double * auxc = malloc(sizeof(double)*4);
    double * aux = malloc(sizeof(double)*4);

    mi=0;                       //T direction
    setquadv(&ami[0],mi);

    for(t=0 ; t<Nt ; t++){
        for(x=0 ; x<N ; x++){
            for(y=0 ; y<N ;y++){
                for(z=0 ; z<N ; z++){
                    for(ni=1 ; ni<4 ; ni++){
                        setquadv(&ani[0],ni);
                        //for(mi=ni+1 ; mi<4 ; mi++){
                            k=0;
                            j=0;
                            setidv(aux);

                            for(k=0 ; k<T ; k++){
								mmprodv(aux, aux, getU(lattice, (t+k)%Nt , x , y , z ,mi) );
							}
							k--;

                            for(j=0 ; j<R ; j++){
								mmprodv(aux , aux , getU(lattice, (t+(k+1)*ami[0])%Nt , (x+(k+1)*ami[1]+ j*ani[1])%N ,
										(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
							}
							j--;

                            for( ; k > -1 ; k--){
								hermcv( auxc , getU(lattice, (t+k+ (j+1)*ani[0])%Nt , (x+ (j+1)*ani[1])%N ,
									(y+ (j+1)*ani[2])%N , (z+ (j+1)*ani[3])%N , mi ) );
								mmprodv(aux, aux , auxc);
							}
							k++;

							for( ; j > -1 ; j--){
								hermcv( auxc , getU(lattice, t , (x+ j*ani[1])%N ,
									(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
								mmprodv(aux , aux , auxc);
							}

                            sum += tracev(aux);
                        //}
                    }
                }
            }
        }
    }
    free(auxc);
	free(aux);
	return sum/(2*3.0*totalV);
}

double polyakovLoopMeasure(double * lattice){
    double sum = 0;
    int temp;
    int k,t,x,y,z;
    double * aux = malloc(sizeof(double)*4);

    t=0;
    for(x=0 ; x<N ; x++){
        for(y=0 ; y<N ;y++){
            for(z=0 ; z<N ; z++){
                setidv(aux);
                for(k=0 ; k<Nt ; k++){
                    mmprodv(aux, aux, getU(lattice, t+k , x, y, z, 0));
                }

                sum += tracev(aux);
            }
        }
    }
    free(aux);
    return sum/(2.0*spatialV);
}

double measure_ghostp(double * lattice, double *** G, double * kz, int kSize, int j_Ncf){
    clock_t dtime = clock();
    int i;
    double complex * chute 		= malloc(sizeof(double complex)*colorV);
    double complex * source 	= malloc(sizeof(double complex)*colorV);
    double p[4] = {0,0,0,0};

    printf("\nGhost processing initiated.\n");
    for(i=0;i<kSize;i++){
        kz[i] = (double)(i+1)/(N);

        //define source to plane wave
        //can add other colors here
        p[3] = kz[i];                   //set suitable 4vector, considering symmetric matrix
        setplanewave(source,p,0);

        //define initial value for CG ??????? should be zero???????
        setzerovc(chute,colorV);

        //getOrthogonalSpace(source, source, lattice);    //eliminating null vectors
        fpcgmr(chute, lattice , source, colorV);	//chute = (M^-1) . source = (M^-1) . M . planewave

        (*G)[i][j_Ncf] = inprodvc(source,chute,colorV)/(double)(totalV);
        printf("G[%.4lf]=%lf\n",kz[i],(*G)[i][j_Ncf]);
    }
    return(((double)(clock() - dtime))/(CLOCKS_PER_SEC));
}

double measure_gluonp(double * lattice, double *** D, double * kz, int kSize, int j_Ncf){
    //to pass a pointer to a pointer in C I must declare it as a pointer to a pointer to pointer
    double gluonP(double * lattice, double kz){
        int t,x,y,z,a,mi;
        double D = 0e0;

        if(kz>0){
            double costemp, sintemp, fouriercos, fouriersin;
            for(mi=0;mi<4;mi++){
                for(a=1;a<4;a++){
                    fouriercos = 0e0;
                    fouriersin = 0e0;
                    for(z=0;z<N;z++){
                        costemp = cos(2*M_PI*kz*z);
                        sintemp = sin(2*M_PI*kz*z);
                        for(t=0;t<Nt;t++){
                            for(x=0;x<N;x++){
                                for(y=0;y<N;y++){
                                    fouriercos += getU(lattice,t,x,y,z,mi)[a]*costemp;
                                    fouriersin += getU(lattice,t,x,y,z,mi)[a]*sintemp;
                                }
                            }
                        }
                    }
                    D += pow(fouriercos,2) + pow(fouriersin,2);
                }
            }
            return(D/(9.0*N*N*N*N));
        }
        else{
            for(mi=0;mi<4;mi++){
                for(a=1;a<4;a++){
                    double temp=0e0;
                    for(z=0;z<N;z++){
                        for(t=0;t<Nt;t++){
                            for(x=0;x<N;x++){
                                for(y=0;y<N;y++){
                                    temp += getU(lattice,t,x,y,z,mi)[a];
                                }
                            }
                        }
                    }
                    D += pow(temp,2);
                }
            }
            return(D/(12.0*N*N*N*Nt));
        }
    }

    clock_t dtime = clock();
    int i;
    for(i=0;i<kSize;i++){
        kz[i] = (double)i/N;
        (*D)[i][j_Ncf] = gluonP(lattice,kz[i]);        //store same momentum continuously
        printf("\n%.3lf %e\n",kz[i],(*D)[i][j_Ncf]);
    }
    return(((double)(clock() - dtime))/(CLOCKS_PER_SEC));
}

void spatialSmearing3(double * outlattice, double * inlattice, int n_smear, double alpha_smear){
    double updateStaple(double * lattice , int t, int x, int y, int z, int mi , double * _staple){
    	int i,ni;
        int ani[4], ami[4];
    	double * auxp = malloc((sizeof(double))*4);
    	double * aux1 = malloc((sizeof(double))*4);
    	double * aux2 = malloc((sizeof(double))*4);
    	double * aux3 = malloc((sizeof(double))*4);

        setquadv(&ami[0],mi);
    	setzerov( _staple ); 		//now staple should be a zero DxD matrix
        for(ni=1 ; ni<4 ; ni++){
            if(ni != mi){
    			setquadv(&ani[0],ni);

                copyv(aux1 ,  getU( lattice,t , x , y , z , ni) );
    			copyv( aux2 , getU( lattice, (t+ani[0])%Nt , (x+ani[1])%N , (y+ani[2])%N , (z+ani[3])%N , mi) );
    			hermcv( aux3 , getU( lattice, (t+ami[0])%Nt , (x+ami[1])%N , (y+ami[2])%N , (z+ami[3])%N , ni)  );

    			mmprodv(auxp, aux1 , aux2 );
    			mmprodv(auxp, auxp , aux3);

    			sumv(_staple , _staple , auxp);


    			hermcv(aux1 , getU( lattice, getStepT(t,-ani[0]), getStep(x,-ani[1]) , getStep(y,-ani[2]) , getStep(z,-ani[3]) , ni));
    			copyv(aux2 ,getU( lattice, getStepT(t,-ani[0]) , getStep(x,-ani[1]) , getStep(y,-ani[2]) , getStep(z,-ani[3]) , mi) );
    			copyv(aux3 , getU( lattice, getStepT(t,ami[0]-ani[0]), getStep(x,ami[1]-ani[1]) , getStep(y,ami[2]-ani[2]) , getStep(z,ami[3]-ani[3]) , ni));

    			mmprodv(auxp , aux1 , aux2);
    			mmprodv(auxp , auxp , aux3);

    			sumv(_staple , _staple , auxp);
    		}
    	}
    	free(auxp);
    	free(aux1);
    	free(aux2);
    	free(aux3);
    }

    void spatialSmearStep(double * outlattice, double * inlattice, double alpha_smear){
        int t,x,y,z,i,ai[4],j,aj[4];
        double * aux = malloc(sizeof(double)*4);
        double * herm = malloc(sizeof(double)*4);
        double * sum = malloc(sizeof(double)*4);
        for(t=0;t<Nt;t++){
            for(x=0;x<N;x++){
                for(y=0;y<N;y++){
                    for(z=0;z<N;z++){
                        for(i=1;i<4;i++){
                            cmprodv(sum, alpha_smear, getU(inlattice,t,x,y,z,i));
                            updateStaple(inlattice, t,x,y,z,i,aux);
                            //hermcv(aux,aux);
                            sumv(sum,sum,aux);
                            reunitv(sum);
                            copyv(getU(outlattice,t,x,y,z,i) , sum);
                        }
                    }
                }
            }
        }
        free(sum);
        free(herm);
        free(aux);
    }
    int i;
    copyl(outlattice,inlattice);
    double * templattice = malloc(sizeof(double)*dimLattice);
    for(i=0;i<n_smear;i++){
        copyl(templattice,outlattice);
        spatialSmearStep(outlattice,templattice,alpha_smear);
    }
    free(templattice);
}
