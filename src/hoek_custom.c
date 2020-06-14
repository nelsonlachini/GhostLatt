#include "hoek_custom.h"
#include "global.h"

//GLOBAL VARIABLES

double map(double n, int N_ins){
	return(n - N_ins/2e0);
}

double mapR(double y){
	return y;
}

void imposePeriodicity(double * lattice, double * x0, int N_ins){
	int mu,i;
	int xv[4];
	double * aux = malloc(sizeof(double)*4);
	int nmu[4];

	for(nmu[0]=0 ; nmu[0]<N_ins ; nmu[0]++){
		for(nmu[1]=0 ; nmu[1]<N_ins ; nmu[1]++){
			for(nmu[2]=0 ; nmu[2]<N_ins ; nmu[2]++){
				for(nmu[3]=0 ; nmu[3]<N_ins ; nmu[3]++){
					for(mu=0 ; mu<4 ; mu++){
							//if( ((nmu[0]==N_ins-1) && mu==0) || ((nmu[1]==N_ins-1) && mu==1) || ((nmu[2]==N_ins-1) && mu==2) || ((nmu[3]==N_ins-1) && mu==3) ){
							if( ((nmu[0]==N_ins-1)) || ((nmu[1]==N_ins-1)) || ((nmu[2]==N_ins-1)) || ((nmu[3]==N_ins-1)) ){
								for(i=0;i<4;i++){
									if(nmu[i] == N_ins-1){
										xv[i] = 0;
									}
									else{
										xv[i] = nmu[i];
									}
								}
								mmprodv(aux, getU(lattice,nmu[0]+x0[0],nmu[1]+x0[1],nmu[2]+x0[2],nmu[3]+x0[3],mu), getU(lattice,xv[0]+x0[0],xv[1]+x0[1],xv[2]+x0[2],xv[3]+x0[3],mu) );
								su2sroot(getU(lattice,nmu[0]+x0[0],nmu[1]+x0[1],nmu[2]+x0[2],nmu[3]+x0[3],mu) , aux);
						}
					}
				}
			}
		}
	}
	free(aux);
}

double getStepBorder( double coord , double dcoord , int N_ins){
	//return the lattice coordinate given the increment(useful when there is a chance to get a negative coordinate)
	double aux = modulus(coord+dcoord,N_ins);
	if(aux<0)
		return N_ins+aux;
	return aux;
}

double updateStapleBorder(double * lattice , int t, int x, int y, int z, int mi , double * _staple, double * x0, int N_ins){
	int i,ni;
	int ami[4] = {0,0,0,0};
	int ani[4];
	ami[mi] = 1;
	double * auxp = malloc((sizeof(double))*4);
	double * aux1 = malloc((sizeof(double))*4);
	double * aux2 = malloc((sizeof(double))*4);
	double * aux3 = malloc((sizeof(double))*4);

	setzerov( _staple ); 		//now staple should be a zero DxD matrix
	for(ni=0 ; ni<4 ; ni++){
		if(ni != mi){
			for(i=0 ; i<4 ; i++)	ani[i] = 0;
			ani[ni] = 1 ;

			copyv(aux1 , getU( lattice, x0[0]+getStepBorder(t,ami[0],N_ins) , x0[1]+getStepBorder(x,ami[1],N_ins) , x0[2]+getStepBorder(y,ami[2],N_ins) , x0[3]+getStepBorder(z,ami[3],N_ins), ni) );
			hermcv( aux2 , getU( lattice, x0[0]+getStepBorder(t,ani[0],N_ins) , x0[1]+getStepBorder(x,ani[1],N_ins) , x0[2]+getStepBorder(y,ani[2],N_ins) , x0[3]+getStepBorder(z,ani[3],N_ins) , mi) );
			hermcv( aux3 , getU( lattice,t , x , y , z , ni) );

			mmprodv(auxp, aux1 , aux2 );
			mmprodv(auxp, auxp , aux3);

			sumv(_staple , _staple , auxp);

			hermcv(aux1 , getU( lattice, x0[0]+getStepBorder(t,ami[0]-ani[0],N_ins), x0[1]+getStepBorder(x,ami[1]-ani[1],N_ins) , x0[2]+getStepBorder(y,ami[2]-ani[2],N_ins) , x0[3]+getStepBorder(z,ami[3]-ani[3],N_ins) , ni));
			hermcv(aux2 ,getU( lattice, x0[0]+getStepBorder(t,-ani[0],N_ins) , x0[1]+getStepBorder(x,-ani[1],N_ins) , x0[2]+getStepBorder(y,-ani[2],N_ins) , x0[3]+getStepBorder(z,-ani[3],N_ins) , mi) );
			copyv( aux3 , getU( lattice, x0[0]+getStepBorder(t,-ani[0],N_ins) , x0[1]+getStepBorder(x,-ani[1],N_ins) , x0[2]+getStepBorder(y,-ani[2],N_ins) , x0[3]+getStepBorder(z,-ani[3],N_ins) , ni));

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

void coolingInstantonBorderStep(double * lattice, double * x0, int N_ins , int border_length){
	int t,x,y,z,mi;
	double * staple = malloc(sizeof(double)*4);
	double inf = border_length, sup = N-1-border_length;
	for(t=0 ; t<N_ins ; t++){
		for(x=0 ; x<N_ins ; x++){
			for(y=0 ; y<N_ins ; y++){
				for(z=0 ; z<N_ins ; z++){
					for(mi=0 ; mi<4 ; mi++){
						if( ((t<inf) || (t>sup)) || ((x<inf) || (x>sup)) || ((y<inf) || (y>sup)) || ((z<inf) || (z>sup)) ){
							updateStapleBorder(lattice ,t+x0[0],x+x0[1],y+x0[2],z+x0[3], mi, staple, x0, N_ins);
							hermcv(staple , staple);
							cmprodv( getU(lattice,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mi), 1e0/sqrt(detv(staple)) , staple);
						}
					}
				}
			}
		}
	}
	free(staple);
}

double coolInstantonBorder(double * lattice, double * x0, int N_ins, int Ncool, int border_length){
	clock_t stime = clock();
	int i;
	for(i=0;i<Ncool;i++){
		coolingInstantonBorderStep(lattice, x0, N_ins , border_length);
	}
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

void amuregV(double * amu, double t, double x, double y, double z, int mu, double R, int N_ins){
	double xv[4] = {map(t,N_ins),map(x,N_ins),map(y,N_ins),map(z,N_ins)};
	double x2 = xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2] + xv[3]*xv[3];
	double rho = mapR(R);

	amu[0] = 0e0;

	if(mu==0){
		amu[1] = -xv[1];
		amu[2] = -xv[2];
		amu[3] = -xv[3];
	}
	else if(mu==1){
		amu[1] = xv[0];
		amu[2] = xv[3];
		amu[3] = -xv[2];
	}
	else if(mu==2){
		amu[1] = -xv[3];
		amu[2] = xv[0];
		amu[3] = xv[1];
	}
	else if(mu==3){
		amu[1] = xv[2];
		amu[2] = -xv[1];
		amu[3] = xv[0];
	}

	cmprodv(amu,1e0/(x2+rho*rho),amu);
}

void creategV(double * g, double t, double x, double y, double z, int N_ins){
	double xv[4] = {map(t,N_ins),map(x,N_ins),map(y,N_ins),map(z,N_ins)};
	double xnorm = sqrt(xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2] + xv[3]*xv[3]);
	int i;

	for(i=0;i<4;i++){
		g[i] = xv[i];
	}
	cmprodv(g,1e0/xnorm,g);
}

void umuregV(double * umu, double t, double x, double y, double z, int mu, int partition_number,double R, int N_ins){
	int i;
	double anorm;
	double * aux = malloc(sizeof(double)*4);
	double * amu = malloc(sizeof(double)*4);
	double emu[4] = {0,0,0,0};
	double partition_size = 1e0/partition_number;
	setzerov(amu);

	amuregV(aux,t,x,y,z,mu,R,N_ins);
	sumv(amu,amu,aux);

	for(i=0;i<partition_number;i++){
		emu[mu] += partition_size;
		amuregV(aux,t+emu[0],x+emu[1],y+emu[2],z+emu[3],mu,R,N_ins);
		sumv(amu,amu,aux);
	}

	cmprodv(amu,1e0/(partition_number+1),amu);

	anorm = sqrt(amu[1]*amu[1] + amu[2]*amu[2] + amu[3]*amu[3]);

	if(anorm!=0){
		umu[0] = cos(anorm);
		for(i=1;i<4;i++){
			umu[i] = amu[i]*sin(anorm)/anorm;
		}
	}
	else{
		setidv(umu);
	}

	free(aux);
	free(amu);
}

void singtrans(double * uout, double * uin, double t, double x, double y, double z, int mu, int N_ins){
	int emu[4];
	double * g = malloc(sizeof(double)*4);
	double * gmu = malloc(sizeof(double)*4);
	setquadv(&emu[0],mu);

	creategV(g,t,x,y,z,N_ins);

	creategV(gmu,t+emu[0],x+emu[1],y+emu[2],z+emu[3],N_ins);
	hermcv(gmu,gmu);

	mmprodv(uin,g,uin);
	mmprodv(uout,uin,gmu);

	free(g);
	free(gmu);
}

double teperInstanton(double * U, int partition_number, double R, double c0, double * x0, int N_ins, int N_cool){

	clock_t stime = clock();
	int t,x,y,z,mu,i;
	int emu[4];
	//c0=0;
	for(t=0 ; t<N_ins ; t++){
		for(x=0 ; x<N_ins ; x++){
			for(y=0 ; y<N_ins ; y++){
				for(z=0 ; z<N_ins ; z++){
					for(mu=0 ; mu<4 ; mu++){
						//printpos(t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu);
						//getchar();
						umuregV(getU(U,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu),t-c0,x-c0,y-c0,z-c0,mu,partition_number,R,N_ins);
						singtrans(getU(U,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu),getU(U,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu),t-c0,x-c0,y-c0,z-c0,mu,N_ins);
					}
				}
			}
		}
	}

	imposePeriodicity(U , x0, N_ins);
	//coolLattice(U, N_cool);
	for(i=0;i<N_ins/2;i=i+2){
		coolInstantonBorder(U, x0, N_ins, N_cool, i);
	}
}

double teperAntiInstanton(double * U, int partition_number, double R, double c0, double * x0, int N_ins, int N_cool){
	clock_t stime = clock();
	int t,x,y,z,mu,i;
	int emu[4];
	for(t=0 ; t<N_ins ; t++){
		for(x=0 ; x<N_ins ; x++){
			for(y=0 ; y<N_ins ; y++){
				for(z=0 ; z<N_ins ; z++){
					for(mu=0 ; mu<4 ; mu++){
						umuregV(getU(U,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu),t-c0,x-c0,y-c0,z-c0,mu,partition_number,R,N_ins);
						singtrans(getU(U,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu),getU(U,t+x0[0],x+x0[1],y+x0[2],z+x0[3],mu),t-c0,x-c0,y-c0,z-c0,mu,N_ins);
					}
				}
			}
		}
	}

	imposePeriodicity(U , x0, N_ins);
	//coolLattice(U, N_cool);
	for(i=0;i<N_ins/2;i=i+2){
		coolInstantonBorder(U, x0, N_ins, N_cool, i);
	}
}
