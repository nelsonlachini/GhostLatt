#include "thermalization_hb.h"
#include "global.h"

double updateStaple(double * lattice , int t, int x, int y, int z, int mi , double * _staple){
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

			copyv(aux1 , getU( lattice, (t+ami[0])%Nt , (x+ami[1])%N , (y+ami[2])%N , (z+ami[3])%N, ni) );
			hermcv( aux2 , getU( lattice, (t+ani[0])%Nt , (x+ani[1])%N , (y+ani[2])%N , (z+ani[3])%N , mi) );
			hermcv( aux3 , getU( lattice,t , x , y , z , ni) );

			mmprodv(auxp, aux1 , aux2 );
			mmprodv(auxp, auxp , aux3);

			sumv(_staple , _staple , auxp);

			hermcv(aux1 , getU( lattice, getStepT(t,ami[0]-ani[0]), getStep(x,ami[1]-ani[1]) , getStep(y,ami[2]-ani[2]) , getStep(z,ami[3]-ani[3]) , ni));
			hermcv(aux2 ,getU( lattice, getStepT(t,-ani[0]) , getStep(x,-ani[1]) , getStep(y,-ani[2]) , getStep(z,-ani[3]) , mi) );
			copyv( aux3 , getU( lattice, getStepT(t,-ani[0]) , getStep(x,-ani[1]) , getStep(y,-ani[2]) , getStep(z,-ani[3]) , ni));

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

void su2Bath( double * X , double a ){
	int i;
	double lambda2, length, aux, r[4];
	do{
		for(i=1 ; i<4 ; i++){
			r[i] = ran0(global_seed);
		}
		lambda2 = (log(r[3]) + log(r[2])*pow(cos(2*M_PI*r[1]),2) )/(-2*a*beta);
	}while( pow(ran0(global_seed),2) > 1 - lambda2);
	r[0] = 1 - 2*lambda2;

	length = sqrt(1-r[0]*r[0]);

	do{
		for(i=1 ; i<4 ; i++){
			r[i] = 2*ran0(global_seed) - 1;
		}
		aux = r[1]*r[1] + r[2]*r[2] + r[3]*r[3];
	}while( aux > 1 );
	for(i=1 ; i<4 ; i++){
		r[i] = r[i]*length/sqrt(aux);
	}

	for(i=0 ; i < 4 ; i++){
			*getelv(X , i) = r[i];
	}
}

void microCanonicalStep(double * link, double * _staple){
	double * aux = malloc( sizeof(double)*4 );
	double * staple_aux = malloc( sizeof(double)*4 );

	hermcv(staple_aux, _staple);

	mmprodv(aux, link, _staple);
	cmprodv(staple_aux, tracev(aux)/detv(_staple), staple_aux);

	cmprodv(aux, -1e0, link);

	sumv(link, staple_aux, aux);

	free(aux);
	free(staple_aux);
}

void microCanonicalUpdate(double * lattice){
	int t,x,y,z,mi;
	double * staple = malloc( sizeof(double)*4 );
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						updateStaple(lattice,t,x,y,z,mi,staple);
						microCanonicalStep(getU(lattice,t,x,y,z,mi),staple);
					}
				}
			}
		}
	}
	free(staple);
}

void heatBathStep(double * link, double a, double * Vdagger){
	double * X = malloc( sizeof(double)*4 );
	su2Bath(X , a );
	mmprodv(link , X , Vdagger);
	free(X);
}

void heatbathUpdate(double * lattice){
	int t,x,y,z,mi;
	double * staple = malloc( sizeof(double)*4 );
	double * Vdagger= malloc( sizeof(double)*4 );
	double a;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						updateStaple(lattice ,t,x,y,z,mi , staple);
						hermcv(Vdagger , staple);
						a = sqrt(detv(staple));
						cmprodv(Vdagger , 1/a , Vdagger);
						heatBathStep(getU(lattice,t,x,y,z,mi), a, Vdagger);
					}
				}
			}
		}
	}
	free(staple);
	free(Vdagger);
}

void overrelaxationStepSU2(double * lattice , int t, int x, int y, int z, int mi, double * Vdagger){
	double * Udagger = malloc((sizeof(double))*4);
	hermcv( Udagger , getU( lattice,t,x,y,z,mi));
	mmprodv( Udagger , Udagger , Vdagger);
	mmprodv( Udagger , Vdagger , Udagger);
	copyv( getU( lattice,t,x,y,z,mi) ,Udagger);
	free(Udagger);
}

double updateLattice(double * lattice, int N_hb, int N_mc){
	clock_t dtime = clock();
	int i;
	//int N_hb = 1; //default heat bath value

	//MICROCANONICAL UPDATE
	for(i=0;i<N_mc;i++)
		microCanonicalUpdate(lattice);

	//HEAT BATH UPDATE
	for(i=0;i<N_hb;i++)
		heatbathUpdate(lattice);

	return( ((double)dtime)/CLOCKS_PER_SEC );
}

double thermalizeLattice(double * lattice , int n , int n_hb , int n_over){
	int i;
	clock_t dtime = clock();

	printf("\nThermalization initiated (%d sweeps : %d heat-baths/sweep and %d overrelaxations/sweep)\n",n,n_hb,n_over);
	fflush(stdout);
	progress_panel(0,n);
	for(i=0 ; i < n ; i++){
		updateLattice(lattice , n_hb,  n_over);
		progress_panel(i+1,n);
	}
	printf("\nLattice thermalized with %d updates.\n",n);
	fflush(stdout);

	dtime = clock() - dtime;
	return( ((double)dtime)/CLOCKS_PER_SEC );
}
