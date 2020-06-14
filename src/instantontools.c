#include "instantontools.h"
#include "global.h"

#include "thermalization_hb.h"

//////////////////////////tools

void getImagv(double * Mout, double * vin){
	//vin is an SU2 vector (mapped from the complex matrix)
	int dim = 2;
	double v3,v1;
	v3 = *getelv(vin,3);
	v1 = *getelv(vin,1);
	*getelmr(Mout,0,0,dim) = v3;
	*getelmr(Mout,0,1,dim) = v1;
	*getelmr(Mout,1,0,dim) = v1;
	*getelmr(Mout,1,1,dim) =-v3;
}

void plaquetteProd(double * Vout, double * lattice , int t, int x, int y, int z, int mu, int nu){

	int i, anu[4], amu[4];
	setquadv(&amu[0],mu);
	setquadv(&anu[0],nu);
	double xv[4] = {t,x,y,z};

	//if(xv[mu]==N-1 || xv[nu]==N-1){
		setidv(Vout);
	//}
	//else{
		double * aux = malloc(sizeof(double)*4);
		setidv(Vout);
		mmprodv(Vout, Vout, getU(lattice,t,x,y,z,mu));
		//printm(Vout);
		//printm(getU(lattice,t,x,y,z,mu));

		mmprodv(Vout, Vout, getU(lattice,getStep(t,amu[0]),getStep(x,amu[1]),getStep(y,amu[2]),getStep(z,amu[3]),nu));
		//printm(Vout);
		//printm(getU(lattice,getStep(t,amu[0]),getStep(x,amu[1]),getStep(y,amu[2]),getStep(z,amu[3]),nu));

		hermcv(aux, getU(lattice,getStep(t,anu[0]),getStep(x,anu[1]),getStep(y,anu[2]),getStep(z,anu[3]),mu) );
		mmprodv(Vout, Vout, aux);
		//printm(Vout);
		//printm(getU(lattice,getStep(t,anu[0]),getStep(x,anu[1]),getStep(y,anu[2]),getStep(z,anu[3]),mu));

		hermcv(aux, getU(lattice,t,x,y,z,nu) );
		mmprodv(Vout, Vout, aux);
		//printm(Vout);
		//printm(getU(lattice,t,x,y,z,nu));

		free(aux);
	//}
}

void plaquetteProdG(double * Vout, double * lattice , int t, int x, int y, int z, int mu, int nu){
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	double * aux3 = malloc(sizeof(double)*4);
	double * aux4 = malloc(sizeof(double)*4);

	int i, anu[4], amu[4];
	setquadv(&amu[0],mu);
	setquadv(&anu[0],nu);

	setidv(Vout);

	if(mu>0){
		copyv(aux1,getU(lattice,t,x,y,z,mu) );
		hermcv(aux3, getU(lattice,getStep(t,anu[0]),getStep(x,anu[1]),getStep(y,anu[2]),getStep(z,anu[3]) ,mu) );
	}
	else{
		hermcv(aux1, getU(lattice,getStep(t,amu[0]),getStep(x,amu[1]),getStep(y,amu[2]),getStep(z,amu[3]),-mu) );
		copyv(aux3, getU(lattice,getStep(t,amu[0]+anu[0]),getStep(x,amu[1]+anu[1]),getStep(y,amu[2]+anu[2]),getStep(z,amu[3]+anu[3]),-mu) );
	}

	if(nu>0){
		copyv(aux2, getU(lattice,getStep(t,amu[0]),getStep(x,amu[1]),getStep(y,amu[2]),getStep(z,amu[3]),nu) );
		hermcv(aux4, getU(lattice,t,x,y,z,nu) );
	}
	else{
		hermcv(aux2, getU(lattice,getStep(t,amu[0]+anu[0]),getStep(x,amu[1]+anu[1]),getStep(y,amu[2]+anu[2]),getStep(z,amu[3]+anu[3]) ,-nu) );
		copyv(aux4, getU(lattice,getStep(t,anu[0]),getStep(x,anu[1]),getStep(y,anu[2]),getStep(z,anu[3]) ,-nu) );
	}

	mmprodv(Vout, Vout, aux1);

	mmprodv(Vout, Vout, aux2);

	mmprodv(Vout, Vout, aux3);

	mmprodv(Vout, Vout, aux4);

	free(aux1);
	free(aux2);
	free(aux3);
	free(aux4);
}

void cloverProd(double * Vout, double * lattice , int t, int x, int y, int z, int mu, int nu){
	double * aux = malloc(sizeof(double)*4);
	int i, anu[4], amu[4];
	setquadv(&amu[0],mu);
	setquadv(&anu[0],nu);

	setidv(Vout);

	mmprodv(Vout, Vout, getU(lattice,t,x,y,z,mu));
	mmprodv(Vout, Vout, getU(lattice,getStep(t,amu[0]),getStep(x,amu[1]),getStep(y,amu[2]),getStep(z,amu[3]),mu));

	mmprodv(Vout, Vout, getU(lattice,getStep(t,2*amu[0]),getStep(x,2*amu[1]),getStep(y,2*amu[2]),getStep(z,2*amu[3]),nu));
	mmprodv(Vout, Vout, getU(lattice,getStep(t,2*amu[0]+anu[0]),getStep(x,2*amu[1]+anu[1]),getStep(y,2*amu[2]+anu[2]),getStep(z,2*amu[3]+anu[3]),nu));


	hermcv(aux, getU(lattice,getStep(t,amu[0]+2*anu[0]),getStep(x,amu[1]+2*anu[1]),getStep(y,amu[2]+2*anu[2]),getStep(z,amu[3]+2*anu[3]),mu));
	mmprodv(Vout, Vout, aux);
	hermcv(aux, getU(lattice,getStep(t,2*anu[0]),getStep(x,2*anu[1]),getStep(y,2*anu[2]),getStep(z,2*anu[3]),mu));
	mmprodv(Vout, Vout, aux);

	hermcv(aux, getU(lattice,getStep(t,anu[0]),getStep(x,anu[1]),getStep(y,anu[2]),getStep(z,anu[3]),nu));
	mmprodv(Vout, Vout, aux);
	hermcv(aux, getU(lattice,t,x,y,z,nu) );
	mmprodv(Vout, Vout, aux);

	free(aux);
}

double integrate(double * indens){
	double Qtemp = 0e0;
	int t,x,y,z;
	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp += *(indens + t*N*N*N + x*N*N + y*N + z);
				}
			}
		}
	}
	return(Qtemp);
}

void integrate1d(double * outdenst, double * indens){
	double Qtemp;
	int t,x,y,z;
	double * tempdenst = malloc(sizeof(double)*N);

	for(t=0;t<N;t++){
		Qtemp = 0e0;
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp += *(indens + t*N*N*N + x*N*N + y*N +z);
				}
			}
		}
		tempdenst[t] = Qtemp;
	}
	copyvr(outdenst, tempdenst , N);
	free(tempdenst);
}

void integrate2d(double * outdenst, double * indens){
	double Qtemp;
	int t,x,y,z;
	double * tempdenst = malloc(sizeof(double)*N*N);

	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			Qtemp = 0e0;
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp += *(indens + t*N*N*N + x*N*N + y*N + z);
				}
			}
			*(tempdenst + t*N + x) = Qtemp;
		}

	}
	copyvr(outdenst, tempdenst , N*N);
	free(tempdenst);
}

///////////////////////////smoothings
void coolingStep(double * lattice){
	int t,x,y,z,mi;
	double * staple = malloc(sizeof(double)*4);
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						updateStaple(lattice ,t,x,y,z,mi , staple);
						hermcv( staple , staple);
						cmprodv( getU(lattice,t,x,y,z,mi), 1e0/sqrt(detv(staple)) , staple);
					}
				}
			}
		}
	}
	free(staple);
}

double coolLattice(double * lattice, int Ncool){
	clock_t stime = clock();
	int i;
	for(i=0;i<Ncool;i++){
		coolingStep(lattice);
	}
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

void coolingStepBorder(double * lattice){
	int t,x,y,z,mi;
	double * staple = malloc(sizeof(double)*4);
	double inf = 2, sup = N-2;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						if( ((t<inf) || (t>sup)) || ((x<inf) || (x>sup)) || ((y<inf) || (y>sup)) || ((z<inf) || (z>sup)) ){
							//printpos(t,x,y,z,mi);
							updateStaple(lattice ,t,x,y,z,mi , staple);
							hermcv( staple , staple);
							cmprodv( getU(lattice,t,x,y,z,mi), 1e0/sqrt(detv(staple)) , staple);
						}
					}
				}
			}
		}
	}
	free(staple);
}

double coolLatticeBorders(double * lattice, int Ncool){
	clock_t stime = clock();
	int i;
	for(i=0;i<Ncool;i++){
		coolingStepBorder(lattice);
	}
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

void APEsmearStep(double * lattice, double alphaape){
	int t,x,y,z,mi;
	double * staple = malloc(sizeof(double)*4);
	double * aux = malloc(sizeof(double)*4);
	double * temp = malloc(sizeof(double)*N*N*N*N*4*4);
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						cmprodv(aux,1.0-alphaape,getU(lattice,t,x,y,z,mi));

						updateStaple(lattice, t, x, y, z, mi, staple);
						cmprodv(staple,alphaape/6e0,staple);

						sumv(getU(temp,t,x,y,z,mi), aux, staple);
						cmprodv(getU(temp,t,x,y,z,mi), 1e0/sqrt(detv(getU(temp,t,x,y,z,mi))) , getU(temp,t,x,y,z,mi));
					}
				}
			}
		}
	}
	copyl(lattice,temp);
	free(staple);
	free(aux);
	free(temp);
}

double APEsmearLattice(double * lattice, double alphaape, int Ncool){
	clock_t stime = clock();
	int i;
	for(i=0;i<Ncool;i++){
		APEsmearStep(lattice,alphaape);
	}
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

/////////////////////action definitions
double calcSwilson(double * lattice, double * sdens){
	clock_t stime = clock();
	double Stemp = 0e0;
	int mu,nu;
	int t,x,y,z;
	double * aux1;
	aux1 = malloc(sizeof(double)*4);

	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Stemp = 0e0;
					for(mu=0;mu<4;mu++){
						for(nu=0;nu<4;nu++){
							if(mu>nu){
								plaquetteProd(aux1,lattice,t,x,y,z,mu,nu);			//aux used as SU2-4vectors
								Stemp += 1e0 - 5e-1*tracev(aux1);
							}
						}
					}
					*(sdens+t*N*N*N + x*N*N + y*N + z) = Stemp;///(2*M_PI*M_PI);
				}
			}
		}
	}
	free(aux1);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

/////////////////////topological charge field-theoretical definitions
double calcQnaive(double * lattice, double * qdens){
	clock_t stime = clock();
	double Qtemp;
	int m,n,r,s;
	int t,x,y,z;
	double * aux1, * aux2;
	aux1 = malloc(sizeof(double)*4);
	aux2 = malloc(sizeof(double)*4);
	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp = 0e0;
					for(m=0;m<4;m++){
						for(n=0;n<4;n++){
							for(r=0;r<4;r++){
								for(s=0;s<4;s++){
									if( (m!=n)&&(m!=r)&&(m!=s)&&(n!=r)&&(n!=s)&&(r!=s)  ){	//all indices are differentn from each other
										plaquetteProd(aux1,lattice,t,x,y,z,m,n);			//aux used as SU2-4vectors
										plaquetteProd(aux2,lattice,t,x,y,z,r,s);

										mmprodv(aux1, aux1, aux2);

										Qtemp += tracev(aux1) * eps4(m,n,r,s);
									}
								}
							}
						}
					}
					*(qdens+t*N*N*N + x*N*N + y*N +z) = -Qtemp/(32*M_PI*M_PI);
				}
			}
		}
	}
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double calcQnaiveimag(double * lattice , double * qdens){	//from Review
	clock_t stime = clock();
	double Qtemp;
	int m,n,r,s;
	int t,x,y,z;
	double * aux1, * aux2;
	aux1 = malloc(sizeof(double)*4);
	aux2 = malloc(sizeof(double)*4);

	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp = 0e0;
					for(m=0;m<4;m++){
						for(n=0;n<4;n++){
							for(r=0;r<4;r++){
								for(s=0;s<4;s++){
									if( (m!=n)&&(m!=r)&&(m!=s)&&(n!=r)&&(n!=s)&&(r!=s)  ){	//all indices are differentn from each other
										plaquetteProd(aux1,lattice,t,x,y,z,m,n);			//aux used as SU2-4vectors
										plaquetteProd(aux2,lattice,t,x,y,z,r,s);

										getImagv(aux1,aux1);								//aux used as double 2x2 matrix
										getImagv(aux2,aux2);

										mmprodr(aux1, aux1, aux2, 2);

										Qtemp += tracer(aux1,2) * eps4(m,n,r,s);
									}
								}
							}
						}
					}
					*(qdens+t*N*N*N + x*N*N + y*N +z) = Qtemp/(32*M_PI*M_PI);
				}
			}
		}
	}
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

//REPAIR DUE TO P.B.C TO BE MADE
double calcQcloverimag(double * lattice , double * qdens){	//from Review
	clock_t stime = clock();
	double Qtemp;
	int m,n,r,s;
	int t,x,y,z;
	double * aux1, * aux2;
	aux1 = malloc(sizeof(double)*4);
	aux2 = malloc(sizeof(double)*4);

	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp = 0e0;
					for(m=0;m<4;m++){
						for(n=0;n<4;n++){
							for(r=0;r<4;r++){
								for(s=0;s<4;s++){
									if( (m!=n)&&(m!=r)&&(m!=s)&&(n!=r)&&(n!=s)&&(r!=s)  ){	//all indices are differentn from each other
										cloverProd(aux1,lattice,t,x,y,z,m,n);				//aux used as SU2-4vectors
										cloverProd(aux2,lattice,t,x,y,z,r,s);

										getImagv(aux1,aux1);								//aux used as double 2x2 matrix
										getImagv(aux2,aux2);

										mmprodr(aux1, aux1, aux2, 2);

										Qtemp += tracer(aux1,2) * eps4(m,n,r,s);
									}
								}
							}
						}
					}
					*(qdens+t*N*N*N + x*N*N + y*N +z) = Qtemp/(16*32*M_PI*M_PI);
				}
			}
		}
	}
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double calcQnaivesym(double * lattice,double * qdens){	//from Vecchia's paper and Rothe
	clock_t stime = clock();
	double Qtemp;
	int m,n,r,s;
	int t,x,y,z;
	double * aux1, * aux2;
	aux1 = malloc(sizeof(double)*4);
	aux2 = malloc(sizeof(double)*4);

	for(t=0;t<N;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					Qtemp = 0e0;
					for(m=-4;m<5;m++){
						for(n=-4;n<5;n++){
							for(r=-4;r<5;r++){
								for(s=-4;s<5;s++){
									if( ((abs(m)!=abs(n))&&(abs(m)!=abs(r))&&(abs(m)!=abs(s))&&(abs(n)!=abs(r))&&(abs(n)!=abs(s))&&(abs(r)!=abs(s))) && (m*n*r*s != 0)  ){	//all indices are differentn from each other and none is zero
										plaquetteProdG(aux1,lattice,t,x,y,z,m,n);
										plaquetteProdG(aux2,lattice,t,x,y,z,r,s);

										mmprodv(aux1, aux1, aux2);

										Qtemp += tracev(aux1) * eps4tilde(m,n,r,s);
									}
								}
							}
						}
					}
					*(qdens+t*N*N*N + x*N*N + y*N +z) = -Qtemp/(16*32*M_PI*M_PI);
				}
			}
		}
	}
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

//AUXILIAR

void saveChargeProfiles(double * Qdensity, char file_name[50], FILE * f){
	int i,j;
	clock_t dtime = clock();
	double * Q1d = malloc(sizeof(double)*N);
	double * Q2d = malloc(sizeof(double)*N*N);
	char aux[20], temp_name[50];
	integrate1d(Q1d,Qdensity);
	for(i=0;i<N;i++)
		printf("\n%lf",Q1d[i]);

	strcpy(temp_name, file_name);
	sprintf(aux , "N%d1d.dat" , N);
	strcat(temp_name , aux);

	f = fopen( temp_name , "w");
	fprintf(f, "N=%d\n" , N);
	fprintf(f, "Qtop=%e\n",integrate(Qdensity));
	fprintf(f,"TOTAL EXECUTION TIME: %lf hours\n",(double)(clock()-dtime)/CLOCKS_PER_SEC*3600.0);
	for(i=0;i<N;i++)
		fprintf(f,"%d	%e\n",i,Q1d[i]);
	fprintf(f,"%d	%e\n",i,Q1d[0]);
	fclose(f);

	integrate2d(Q2d,Qdensity);

	strcpy(temp_name, file_name);
	sprintf(aux , "N%d2d.dat" , N);
	strcat(temp_name , aux);

	f = fopen( temp_name , "w");
	fprintf(f, "N = %d\n" , N);
	fprintf(f, "Qtop=%e\n",integrate(Qdensity));
	fprintf(f,"TOTAL EXECUTION TIME: %lf hours\n",(double)(clock()-dtime)/CLOCKS_PER_SEC*3600.0);
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fprintf(f,"%d	%d	%e\n",i,j,Q2d[i*N+j]);
		}
	}
	fclose(f);

	free(Q1d);
	free(Q2d);
}

double continuum1dChargeDensity(double rho, double x, double x0){
	x0 = (N+x0)/2.0;
	return( (3.0*pow(rho,4.0))/( 4.0*pow((x-x0)*(x-x0) + rho*rho,2.5) ) );
}


/*

void calcTopCtildeslicedGen(double * lattice, double * Qtilde, int slice_coord){
	double Qtemp;
	int m,n,r,s;
	int i,j;
	int xv[4],aux[3];
	double * aux1, * aux2;
	aux1 = malloc(sizeof(double)*4);
	aux2 = malloc(sizeof(double)*4);


	j=0;
	for(i=0;i<4;i++){
		if(i != slice_coord){
			aux[j] = i;
			j++;
		}
	}

	for(xv[slice_coord]=0;xv[slice_coord]<N;xv[slice_coord]++){
		Qtemp = 0E0;

		for(xv[aux[0]]=0;xv[aux[0]]<N;xv[aux[0]]++){
			for(xv[aux[1]]=0;xv[aux[1]]<N;xv[aux[1]]++){
				for(xv[aux[2]]=0;xv[aux[2]]<N;xv[aux[2]]++){
					for(m=0;m<4;m++){
						for(n=0;n<4;n++){
							for(r=0;r<4;r++){
								for(s=0;s<4;s++){

									if( (m!=n)&&(m!=r)&&(m!=s)&&(n!=r)&&(n!=s)&&(r!=s)  ){	//all indices are differentn from each other

										plaquetteProd(aux1,lattice,xv[0],xv[1],xv[2],xv[3],m,n);			//aux used as SU2-4vectors
										plaquetteProd(aux2,lattice,xv[0],xv[1],xv[2],xv[3],r,s);

										mmprodv(aux1, aux1, aux2);

										Qtemp += tracev(aux1) * eps4(m,n,r,s);
									}
								}
							}
						}
					}
				}
			}
		}
		Qtilde[xv[slice_coord]] = Qtemp/(32*M_PI*M_PI);
	}
	free(aux1);
	free(aux2);
}
*/
