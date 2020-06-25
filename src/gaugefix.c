#include "gaugefix.h"
#include "utilities.h"
#include "algebra.h"

//////////////////////////auxiliars
void initg(double * g){
	int t,x,y,z;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setidv(getg(g,t,x,y,z));
				}
			}
		}
	}
}

void reunitg(double * g){
	int t,x,y,z;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					reunitv(getg(g,t,x,y,z));
				}
			}
		}
	}
}

double calcEps(LatticeSU2 * lattice , double * g){
	int t,x,y,z,mi,i;
	double sum = 0.E0;
	double * auxU = malloc(sizeof(double)*4);
	int ami[4];
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			for(i=0 ; i<4 ; i++)	ami[i] = 0;
			ami[mi] = 1;

			hermcv(auxU , getg(g, getStepT(t,+ami[0]) , getStep(x,+ami[1]) , getStep(y,+ami[2]) , getStep(z,+ami[3])) );
			mmprodv(auxU , getU(lattice,t,x,y,z,mi) , auxU);
			mmprodv(auxU , getg(g,t,x,y,z) , auxU);
			sum += tracev(auxU);
		}
	}}}}
	sum /= (2.0*4*totalV);
	return 1-sum;
}

double applyFix(LatticeSU2 * latticeOut, LatticeSU2 * latticeIn , double * g){
	clock_t stime = clock();
	double * auxU = malloc(sizeof(double)*4);
	int t,x,y,z,mi,i,ami[4];
	for(t=0 ; t<latticeIn->Nt ; t++){
	for(x=0 ; x<latticeIn->N ; x++){
	for(y=0 ; y<latticeIn->N ; y++){
	for(z=0 ; z<latticeIn->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			setquadv(&ami[0],mi);
			hermcv(auxU , getg(g, getStepT(t,+ami[0]) , getStep(x,+ami[1]) , getStep(y,+ami[2]) , getStep(z,+ami[3])) );
			mmprodv(auxU , getU(latticeIn,t,x,y,z,mi) , auxU);
			mmprodv(getU(latticeOut,t,x,y,z,mi) , getg(g,t,x,y,z) , auxU);
		}
	}}}}
	free(auxU);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

void copyg(double * gout, double * gin){
	copyvr(gout,gin,totalV*4);
}

double su2distance(LatticeSU2 * lattice, LatticeSU2 * W){
	int t,x,y,z,mi;
	double d = 0e0;
	double * sum = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);

	setzerov(sum);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0;mi<4;mi++){
			cmprodv(aux1,-1,getU(W,t,x,y,z,mi));
			sumv(aux1, getU(lattice,t,x,y,z,mi) , aux1);
			hermcv(aux2,aux1);
			mmprodv(aux1,aux2,aux1);
			sumv(sum,sum,aux1);
		}
	}}}}
	d = tracev(sum)/(8e0*totalV);
	free(sum);
	free(aux1);
	free(aux2);
	return(d);
}

//SERIAL VERSIONS

double los_alamos(LatticeSU2 * lattice , double * g , double * e2){
	int t,x,y,z,mi,i,ami[4];
	double auxN,auxT;
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setquadv(&ami[0],mi);

			hermcv( aux1 , getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getU(lattice,t,x,y,z,mi) , aux2);
			sumv(h , h , aux1);
		}
		auxN = sqrt(detv(h));						//auxN <- sqrt(|h|) = N
		cmprodv(h , 1e0/auxN , h);					//h <- h/N = ~h
		mmprodv( aux1 , getg(g,t,x,y,z) , h);		//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);						//aux1 <- tr(~w) = T
		*e2 += auxN*auxN*(1e0 - auxT*auxT/4e0);

		hermcv(getg(g,t,x,y,z) , h);				//Los Alamos update: g -> ~h^dagger
	}}}}

	*e2 /= (double)(totalV);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double cornell(LatticeSU2 * lattice , double * g, double * e2 , double alpha){
	int t,x,y,z,mi,i;
	int ami[4];
	double auxN,auxT;
	double * id = malloc(sizeof(double)*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0e0;
	setidv(id);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setquadv(&ami[0],mi);

			hermcv( aux1 , getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getU(lattice,t,x,y,z,mi) , aux2);
			sumv(h , h , aux1);
		}
		//reunitarizing expression
		auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
		cmprodv(h , 1.E0/auxN , h);						//h <- h/N = ~h
		mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);							//auxT <- tr(~w) = T
		hermcv(aux1 , aux1);							//aux1 <- ~w^dag

		cmprodv(aux1 , alpha*auxN , aux1);
		cmprodv(aux2 , 1.E0 - alpha*auxN*auxT/2E0 , id);
		sumv(aux1 , aux1 , aux2);

		reunitv(aux1);				//aux1 <- R_update;

		mmprodv(getg(g,t,x,y,z) , aux1 , getg(g,t,x,y,z));	//overrelaxation update


		//analitic expression
		/*
		auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
		cmprodv(h , 1e0/auxN , h);						//h <- h/N = ~h
		mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);							//auxT <- tr(~w) = T

		hermcv(h,h);
		cmprodv(h,alpha*auxN,h);
		cmprodv(aux1, (1-alpha*auxN*auxT/2e0) , getg(g,t,x,y,z) );
		sumv(aux1,aux1,h);

		cmprodv( getg(g,t,x,y,z) , 1e0/sqrt( 1+alpha*alpha*auxN*auxN*(1-auxT*auxT/4e0) ) , aux1);
		*/

		*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
	}}}}

	*e2 /= (double)(totalV);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double overrelaxation(LatticeSU2 * lattice , double * g, double * e2 , double omega){
	int t,x,y,z,mi,i;
	int ami[4];
	double auxN,auxT;
	double * id = malloc(sizeof(double)*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	setidv(id);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setquadv(&ami[0],mi);

			hermcv( aux1 , getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getU(lattice,t,x,y,z,mi) , aux2);
			sumv(h , h , aux1);
		}

		//reunitarizing version
		auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
		cmprodv(h , 1e0/auxN , h);						//h <- h/N = ~h
		mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);							//auxT <- tr(~w) = T

		hermcv(aux1 , aux1);							//aux1 <- ~w^dag

		cmprodv(aux1 , omega , aux1);
		cmprodv(aux2 , (1e0-omega) , id);
		sumv(aux1 , aux1 , aux2);

		reunitv(aux1);
		mmprodv(getg(g,t,x,y,z) , aux1 , getg(g,t,x,y,z));	//overrelaxation update

		//analitc version
		/*
		auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
		cmprodv(h , 1e0/auxN , h);						//h <- h/N = ~h
		mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);							//auxT <- tr(~w) = T
		hermcv(h,h);

		cmprodv(h,omega,h);
		cmprodv(aux1 , 1e0 - omega, getg(g,t,x,y,z) );
		sumv(aux1,h,aux1);

		cmprodv(getg(g,t,x,y,z) , 1e0/sqrt(1+omega*(omega-1)*(2e0-auxT) ) , aux1);
		*/

		*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
	}}}}
	*e2 /= (double)(totalV);

	free(id);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC );
}

double stochastic(LatticeSU2 * lattice , double * g , double * e2, double p){
	int t,x,y,z,mi,i,ami[4];
	double auxN,auxT;
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0e0;
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setquadv(&ami[0],mi);

			hermcv( aux1 , getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getU(lattice,t,x,y,z,mi) , aux2);
			sumv(h , h , aux1);
		}

		//printf("-------------\n h:\n");
		//printm(h);

		auxN = sqrt(detv(h));						//auxN <- sqrt(|h|) = N
		cmprodv(h , 1e0/auxN , h);					//h <- h/N = ~h
		mmprodv( aux1 , getg(g,t,x,y,z) , h);		//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);						//aux1 <- tr(~w) = T
		//printf("\naux1:\n");
		//printm(aux1);

		if( ran0(global_seed) < p ){
			/*
			hermcv(aux1,aux1);
			mmprodv(aux1 , aux1 , aux1);
			mmprodv(getg(g,t,x,y,z) , aux1 , getg(g,t,x,y,z) );
			*/
			hermcv(aux1,h);
			cmprodv(aux1,auxT,aux1);
			cmprodv(aux2,-1,getg(g,t,x,y,z));
			sumv(getg(g,t,x,y,z),aux1,aux2);
		}
		else{
			hermcv(getg(g,t,x,y,z) , h);					//Los Alamos update
		}

		*e2 += auxN*auxN*(1e0 - auxT*auxT/4e0);

		if(isnan(*e2)){
			printf("\n inside stoch e2 = %.2e ; auxN = %.3e ; auxT = %.3e\n",*e2,auxN,auxT);
			getchar();
		}
	}}}}

	*e2 /= (double)(totalV);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

/*
//in development
double fourier(double * U , double * g , double * e2, double alpha){
	void computew(double * wout, double * U, double * g,int t, int x, int y, int z){
		//STOPPED HERE
		int i,mi,ami[4];
		double * h = malloc(sizeof(double)*4);
		double * aux1 = malloc(sizeof(double)*4);
		double * aux2 = malloc(sizeof(double)*4);
		setzerov(h);
		for(mi=0;mi<4;mi++){
			setquadv(&ami[0],mi);

			hermcv( aux1 , getU(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getU(U,t,x,y,z,mi) , aux2);
			sumv(h , h , aux1);
		}
		mmprodv(wout,getg(g,t,x,y,z),h);
		free(h);
		free(aux1);
		free(aux2);
	}

	int t,x,y,z,mi,i,ami[4];
	double auxN,auxT;
	double * w = malloc(sizeof(double)*4);
	double * wforw = malloc(sizeof(double)*4);
	double * wbackw = malloc(sizeof(double)*4);
	double * lap = malloc(sizeof(double)*4);
	double * auxg = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setzerov(lap);
					//compute lattice laplacian of w
					computew(w,U,g,t,x,y,z);
					auxN = sqrt(detv(w));
					auxT = tracev(w)/auxN;
					cmprodv(w,-8e0,w);
					sumv(lap,lap,w);
					for(mi=0;mi<4;mi++){
						setquadv(&ami[0],mi);
						computew(wbackw,U,g,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]));
						computew(wforw,U,g,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]));
						sumv(lap,lap,wbackw);
						sumv(lap,lap,wforw);
					}
					//invert(learn FFT!!!!!!!!)and .(-1) the laplacian components
					for(i=1;i<4;i++){
						//lap[i] = -alpha/lap[i];
					}


					auxg[0] = 1e0;
					for(i=1;i<4;i++){
						auxg[i] = -lap[i];
					}
					//cmprodv

					cmprodv(auxg , 1e0/sqrt( 1e0+alpha*alpha*(lap[1]*lap[1]+lap[2]*lap[2]+lap[3]*lap[3]) ) , auxg);
					mmprodv(getg(g,t,x,y,z) , auxg , getg(g,t,x,y,z));

					*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
				}
			}
		}
	}
	*e2 /= (double)(totalV);
	free(w);
	free(wforw);
	free(wbackw );
	free(auxg);
	free(lap);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}
*/

//AUTOMATIZATION
double fixLatticeStoch(LatticeSU2 * lout, LatticeSU2 * lin, double * g, double p_stoch, double e2tol){
	clock_t stime = clock();
	printf("\nSO Landau Gauge Fixing initiated. Expected tolerance e2tol = %.0e\n",e2tol);
	//double * g 			= malloc(sizeof(double)*totalV*4);
	//initg(g);
	int cont = 0;
	double e2;
	printf("\nStarting: e2 = %.0e | minizing E=%.0e",e2,calcEps(lin , g));
	printf("\niter=%d | e2 = %.0e | minizing E=%lf",cont,e2,calcEps(lin , g));
	do{
		stochastic(lin,g,&e2,p_stoch);
		printf("\riter=%d | e2 = %.0e | minizing E=%lf",cont,e2,calcEps(lin , g));
		cont++;
	}while(e2 > e2tol);
	applyFix(lout,lin,g);
	//compareLattice(lout,lin);
	//getchar();
	printf("\nGauge fixing completed.\n");
	//free(g);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double fixLatticeLosalamos(LatticeSU2 * lout, LatticeSU2 * lin, double * g , double e2tol){
	clock_t stime = clock();
	printf("\n Los Alamos gauge fixing initiated. Expected tolerance e2tol = %.0e\n",e2tol);
	//initg(g);
	int cont = 0;
	double e2;
	printf("\n iter=%d | e2 = %.0e",cont,e2);
	do{
		los_alamos(lin,g,&e2);
		printf("\r iter=%d | e2 = %.0e",cont,e2);
		cont++;
	}while(e2 > e2tol);
	applyFix(lout,lin,g);
	printf("\nGauge fixing completed in %d Los Alamos steps.\n",cont);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double fixLatticeOver(LatticeSU2 * lout, LatticeSU2 * lin, double * g, double omega, double e2tol){
	clock_t stime = clock();
	printf("\nSO Landau Gauge Fixing initiated. Expected tolerance e2tol = %.0e\n",e2tol);
	//double * g 			= malloc(sizeof(double)*totalV*4);
	//initg(g);
	int cont = 0;
	double e2;
	do{
		overrelaxation(lin, g, &e2, omega);
		printf("\riter=%d | e2 = %.0e",cont,e2);
		cont++;
	}while(e2 > e2tol);
	applyFix(lout,lin,g);
	//compareLattice(lout,lin);
	//getchar();
	printf("\nGauge fixing completed.\n");
	//free(g);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

//AUTO CALIBRATION

double calce2(LatticeSU2 * lattice){
	int t,x,y,z,mi,i;
	double sum = 0;
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	double e2 = 0e0,auxN,auxT;
	int ami[4];
	double * g = malloc(sizeof(double)*totalV*4);
	initg(g);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setquadv(&ami[0],mi);

			hermcv(aux1, getU(lattice,t,x,y,z,mi));
			sumv(aux1, aux1, getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi));
			sumv(h,h,aux1);

			hermcv(aux2, getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi));
			sumv(aux2,aux2,getU(lattice,t,x,y,z,mi));
			cmprodv(aux2, -1, aux2);
			sumv(h,h,aux2);
			//hermcv( aux1 , getU(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			//hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			//mmprodv(aux1 , aux1 , aux2);
			//sumv(h , h , aux1);

			//hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			//mmprodv(aux1 , getU(U,t,x,y,z,mi) , aux2);
			//sumv(h , h , aux1);
		}

		//reunitarizing version
		//auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
		//cmprodv(h , 1e0/auxN , h);						//h <- h/N = ~h
		//mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
		//auxT = tracev(aux1);							//auxT <- tr(~w) = T
		hermcv(aux1, h);
		mmprodv(aux1, aux1,h);
		e2 += tracev(aux1);
	}}}}
	free(g);
	free(h);
	free(aux1);
	free(aux2);
	return e2/(totalV);
}

double calibrate_over(LatticeSU2 * lattice , double tol){
	// the parameter omega from overrelaxation is supposed to obey 1<omega<2
	int i;
	double e2;
	double * g = malloc(sizeof(double)*totalV*4);
	double step_size = 0.1;
	double best_omega = 1.05E0;
	double omega_temp = 1.05E0;
	double best_path = 1e20;
	while(omega_temp < 2.0){
		initg(g);
		i=0;
		do{
			overrelaxation(lattice,g,&e2,omega_temp);
			i++;
			printf("\romega=%.2lf --> e2 = %.0e",omega_temp,e2);
		}while(e2 > tol);
		if(i < best_path){
			best_omega = omega_temp;
			best_path = i;
		}
		printf("\romega = %.2lf -> %d steps.\n",omega_temp,i);
		//printf("\n%.2lf	%d\n",omega_temp,i);
		omega_temp += step_size;
	}
	free(g);
	return(best_omega);
}

double calibrate_cornell(LatticeSU2 * lattice , double tol, double step_size, double lower_alpha, double upper_alpha){
	//the alpha parameter from the Cornell method is only supposed to be small and positive
	//so its calibration is more difficult
	int i,j;
	int n_steps = (int)((upper_alpha-lower_alpha)/step_size);
	double e2;
	double best_alpha;
	double temp_alpha;
	double best_path = 1e20;
	double * g = malloc(sizeof(double)*totalV*4);
	printf("\nCornell calibration initiated:\n");
	for(j=0 ; j < n_steps ; j++){
		initg(g);
		temp_alpha = lower_alpha + j*step_size;
		i=0;
		do{
			cornell(lattice,g,&e2,temp_alpha);
			i++;
			printf("\r alpha=%.2lf --> e2 = %.0e",temp_alpha,e2);
		}while(e2 > tol);
		if(i < best_path){
			best_alpha = temp_alpha;
			best_path = i;
		}
		//temp_alpha -= 0.05E0;

		printf(" ---> %d steps.\n",i);
	}
	free(g);
	return(best_alpha);
}

double calibrate_stoc(LatticeSU2 * lattice , double tol){
	// the p parameter of stochastic overrelaxation is supposed to obey 0<p<1
	int i;
	double e2;
	double best_p = 0.85;
	double temp_p = best_p;
	double best_path = 1e20;
	double step_size = 0.075;
	double * g = malloc(sizeof(double)*totalV*4);
	LatticeSU2 * latticetemp = newLatticeSU2(lattice->N,lattice->Nt);
	printf("\nStochastic Overrelaxation calibration initiated:\n");
	while(temp_p > 0e0 ){
		initg(g);
		i=0;
		do{
			stochastic(lattice,g,&e2,temp_p);
			applyFix(latticetemp,lattice,g);
			e2 = calce2(latticetemp);
			i++;
			printf("\r p =%.2lf --> e2 = %.0e",temp_p,e2);
		}while(e2 > tol);
		if(i < best_path){
			best_p = temp_p;
			best_path = i;
		}
		//printf("\n%.2lf	%d\n",temp_p,i);
		printf(" ---> %d steps.\n",i);
		temp_p -= step_size;
	}
	printf("\n Calibration finished with p_stoch=%.2lf",best_p);
	free(g);
	free(latticetemp->U);
	return(best_p);
}


//####################################################IN TEST
/*
double calcFunctional(LatticeSU2 * U){
	int t,x,y,z,mi,i;
	double sum = 0e0;
	int ami[4];
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						sum += tracev( getU(U, t, x, y, z, mi) );
						//printf("\n tr u = %lf", tracev(getU(U, t, x, y, z, mi)) );
					}
				}
			}
		}
	}
	sum /= (8.0*totalV);
	return 1-sum;
}

double * getLatticeReal(double * real , int t , int x , int y , int z){
	return( real + t*N*N*N + x*N*N + y*N + z);
}

double latticeAvgReal(double * var){
	double sum = 0;
	int t,x,y,z;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
				sum += *getLatticeReal(var , t,x,y,z);
				}
			}
		}
	}
	return(sum/pow(N,4));
}*/

/*void computew(double * wout, LatticeSU2 * lattice, double * g,int t, int x, int y, int z){
	//STOPPED HERE
	int i,mi,ami[4];
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	setzerov(h);
	for(mi=0;mi<4;mi++){
		setquadv(&ami[0],mi);

		hermcv( aux1 , getU(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
		hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
		mmprodv(aux1 , aux1 , aux2);
		sumv(h , h , aux1);

		hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
		mmprodv(aux1 , getU(lattice,t,x,y,z,mi) , aux2);
		sumv(h , h , aux1);
	}
	mmprodv(wout,getg(g,t,x,y,z),h);
	free(h);
	free(aux1);
	free(aux2);
}*/

/*
//VECTORIZED VERSIONS
//FORGET THIS
double los_alamosv(double * U , double * g , double * e2){
	int t,x,y,z,mi,i,ami[4];
	double auxN,auxT;
	double * gtemp = malloc(sizeof(double)*totalV*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setzerov(h);
					for(mi=0 ; mi<4 ; mi++){
						for(i=0 ; i<4 ; i++){
							ami[i] = 0;
						}
						ami[mi] = 1;

						hermcv( aux1 , getU(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
						hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
						mmprodv(aux1 , aux1 , aux2);
						sumv(h , h , aux1);

						hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
						mmprodv(aux1 , getU(U,t,x,y,z,mi) , aux2);
						sumv(h , h , aux1);
					}
					auxN = sqrt(detv(h));						//auxN <- sqrt(|h|) = N
					cmprodv(h , 1.E0/auxN , h);					//h <- h/N = ~h
					mmprodv( aux1 , getg(g,t,x,y,z) , h);		//aux1 <- g . ~h = ~w
					auxT = tracev(aux1);						//aux1 <- tr(~w) = T

					hermcv(getg(gtemp,t,x,y,z) , h);				//Los Alamos update to a temporary g

					*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
				}
			}
		}
	}
	*e2 /= (double)(totalV);		//this is the value from the previous iteration
	copyg(g,gtemp);

	free(h);
	free(aux1);
	free(aux2);
	free(gtemp);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double cornellv(double * U , double * g, double * e2 , double alpha){
	int t,x,y,z,mi,i;
	int ami[4];
	double auxN,auxT;
	double * gtemp = malloc(sizeof(double)*totalV*4);
	double * id = malloc(sizeof(double)*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	setidv(id);
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setzerov(h);
					for(mi=0 ; mi<4 ; mi++){
						for(i=0 ; i<4 ; i++){
							ami[i] = 0;
						}
						ami[mi] = 1;

						hermcv( aux1 , getU(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
						hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
						mmprodv(aux1 , aux1 , aux2);
						sumv(h , h , aux1);

						hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
						mmprodv(aux1 , getU(U,t,x,y,z,mi) , aux2);
						sumv(h , h , aux1);
					}
					auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
					cmprodv(h , 1.E0/auxN , h);						//h <- h/N = ~h
					mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
					auxT = tracev(aux1);							//auxT <- tr(~w) = T
					hermcv(aux1 , aux1);							//aux1 <- ~w^dag

					cmprodv(aux1 , alpha*auxN , aux1);
					cmprodv(aux2 , 1.E0 - alpha*auxN*auxT/2E0 , id);
					sumv(aux1 , aux1 , aux2);

					//cmprodv(aux1 , 1.E0/sqrt(1.E0+ omega*(omega-1.E0)*(2.E0-auxT)) , aux1);	// this normalization failed after some time
					reunitv(aux1);				//aux1 <- R_update;

					mmprodv(getg(gtemp,t,x,y,z) , aux1 , getg(g,t,x,y,z));	//overrelaxation update

					*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
				}
			}
		}
	}
	*e2 /= (double)(totalV);
	copyg(g,gtemp);
	free(id);

	free(gtemp);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double overrelaxationv(double * U , double * g, double * e2 , double omega){
	int t,x,y,z,mi,i;
	int ami[4];
	double auxN,auxT;
	double * gtemp = malloc(sizeof(double)*totalV*4);
	double * id = malloc(sizeof(double)*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	setidv(id);
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setzerov(h);
					for(mi=0 ; mi<4 ; mi++){
						for(i=0 ; i<4 ; i++){
							ami[i] = 0;
						}
						ami[mi] = 1;

						hermcv( aux1 , getU(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
						hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
						mmprodv(aux1 , aux1 , aux2);
						sumv(h , h , aux1);

						hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
						mmprodv(aux1 , getU(U,t,x,y,z,mi) , aux2);
						sumv(h , h , aux1);
					}
					auxN = sqrt(detv(h));							//auxN <- sqrt(|h|) = N
					cmprodv(h , 1.E0/auxN , h);						//h <- h/N = ~h
					mmprodv( aux1 , getg(g,t,x,y,z) , h);			//aux1 <- g . ~h = ~w
					auxT = tracev(aux1);							//auxT <- tr(~w) = T
					hermcv(aux1 , aux1);							//aux1 <- ~w^dag

					cmprodv(aux1 , omega , aux1);
					cmprodv(aux2 , (1.E0-omega) , id);
					sumv(aux1 , aux1 , aux2);

					reunitv(aux1);

					mmprodv(getg(gtemp,t,x,y,z) , aux1 , getg(g,t,x,y,z));	//overrelaxation update

					*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
				}
			}
		}
	}
	*e2 /= (double)(totalV);
	copyg(g,gtemp);

	free(gtemp);
	free(id);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC );
}

double stochasticoverv(double * U , double * g , double * e2, double p, long * seed){
	int t,x,y,z,mi,i,ami[4];
	double auxN,auxT;
	double * gtemp = malloc(sizeof(double)*totalV*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	clock_t stime = clock();
	*e2 = 0.E0;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setzerov(h);
					for(mi=0 ; mi<4 ; mi++){
						for(i=0 ; i<4 ; i++){
							ami[i] = 0;
						}
						ami[mi] = 1;

						hermcv( aux1 , getU(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
						hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
						mmprodv(aux1 , aux1 , aux2);
						sumv(h , h , aux1);

						hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
						mmprodv(aux1 , getU(U,t,x,y,z,mi) , aux2);
						sumv(h , h , aux1);
					}

					auxN = sqrt(detv(h));						//auxN <- sqrt(|h|) = N
					cmprodv(h , 1.E0/auxN , h);					//h <- h/N = ~h
					mmprodv( aux1 , getg(g,t,x,y,z) , h);		//aux1 <- g . ~h = ~w
					auxT = tracev(aux1);						//aux1 <- tr(~w) = T

					if( ran0(seed) < p ){
						//hermcv(aux1 , h);								//aux1 <- ~h^dag
						//cmprodv(aux1 , auxT , aux1);					//aux1 <- T . ~h^dag
						//cmprodv(aux2 , -1.E0 , getg(g,t,x,y,z));			//aux2 <- -g
						//sumv(getg(g,t,x,y,z) , aux1 , aux2);			//Stoc update
						hermcv(aux1,aux1);
						mmprodv(aux1 , aux1 , aux1);
						mmprodv(getg(gtemp,t,x,y,z) , aux1 , getg(g,t,x,y,z) );
					}
					else{
						hermcv(getg(gtemp,t,x,y,z) , h);					//Los Alamos update
					}

					*e2 += auxN*auxN*(1.E0 - auxT*auxT/4.E0);
				}
			}
		}
	}
	*e2 /= (double)(totalV);
	copyg(g,gtemp);

	free(gtemp);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}
*/

/*
double calce6(double * lattice , double * g){
	int xni,ni,j;
	double e6 = 0.E0;
	double * aux = malloc(sizeof(double)*4);
	double * U = malloc(sizeof(double)*dimLattice);	//future copy of lattice
	double * Q = malloc(sizeof(double)*(N*4)*4);		//four lists of N matrices
	double * avgQ = malloc(sizeof(double)*4*4);			//four matrices



	//###############################################################-METHODS
	double * getQ(int ni , int xni){
		return(Q + ni*N + xni*4);
	}					//OK: returns the 4-matrix Q[ni,xni]
	double * getAvgQ(int ni){
		return(avgQ + ni*4);
	}						//OK: returns the 4-matrix average <Qni>
	void calcAvgQ(int ni , double * sum){
		int xni;
		setzerov(sum);
		for(xni=0 ; xni<N ; xni++){
			sumv(sum , getQ(ni , xni) , sum);
		}
		cmprodv(sum , 1.E0/(double)N , sum);
	}			//OK: computes the average <Q_ni> and stores it in the input *sum
	void calcQni( int ni  , double * lattice){
		int i,j;
		int x[4], mi[3];
		double * auxd = malloc(sizeof(double)*4);

		//######################################################################-AUXILIAR FUNCTION
		void sumToQni(int mi , int xni){
			hermcv( auxd , getU(lattice, x[0] , x[1] , x[2] , x[3] , mi) );
			cmprodv( auxd , -1.E0 , auxd);
			sumv( auxd , getU(lattice, x[0] , x[1] , x[2] , x[3] , mi) , auxd );
			cmprodv( auxd , 0.5 , auxd );
			sumv( getQ(ni,xni) , getQ(ni,xni) , auxd );
		}
		//########################################################################################

		j=0;
		for(i=0 ; i<4 ; i++){
			if(i != ni){
				mi[j] = i;
				j++;
			}
		}				//OK: this block fills the m[3] with the indices different from ni
		for(x[ni]=0 ; x[ni]<N ; x[ni]++){
			setzerov(getQ(ni,x[ni]));
			for( x[mi[0]]=0 ; x[mi[0]]<N ; x[mi[0]]++ ){
				for( x[mi[1]]=0 ; x[mi[1]]<N ; x[mi[1]]++ ){
					for( x[mi[2]]=0 ; x[mi[2]]<N ; x[mi[2]]++ ){
						sumToQni(mi[0] , x[ni]);
						sumToQni(mi[1] , x[ni]);
						sumToQni(mi[2] , x[ni]);
					}
				}
			}
		}	//APPRENTLY OK: computes the Q_ni(xni) for each xni for the specified ni, and stores it in Q[ni,xni]	]
		free(auxd);
	}		//OK: computes Q_ni(xni) and stores in Q(ni,xni)
	//########################################################################

	applyFix(U , lattice ,g);	//creating a fixed gauge lattice(note that g is the whole transformation): U = g lattice g^dag
	for(ni=0 ; ni<4 ; ni++){
		calcQni(ni, U);
		calcAvgQ(ni, getAvgQ(ni));
	}	//computes Q_ni and the corresponding <Q_ni>

	for(ni=0 ; ni<4 ; ni++){
		//printm(getAvgQ(ni));
		for(j=1 ; j<4 ; j++){
			for(xni=0 ; xni<N ; xni++){
				//printf("Q_%d(%d)_%d = %f ; ",ni,xni,j,*getelv(getQ(ni,xni),j));
				//printf("ratio = %lf" , *getelv(getQ(ni,xni),j)/(*getelv(getAvgQ(ni),j)));
				//printm(getQ(ni,xni));
				cmprodv(aux , -1.E0 , getAvgQ(ni));
				//printm(aux);
				sumv(aux , getQ(ni,xni) , aux);
				//printm(aux);
				//mmprodv(aux, aux , aux);
				//printm(aux);
				//printf("%lf",*getelv(aux,j)/pow(*getelv(getAvgQ(ni),j) , 2) );
				//getchar();
				e6 += pow(*getelv(aux,j),2)/pow(*getelv(getAvgQ(ni),j) , 2) ;
			}
		}
	}	//computes e6

	free(aux);
	free(avgQ);
	free(Q);
	free(U);
	return(e6/(3.0*4*totalV));
}
*/