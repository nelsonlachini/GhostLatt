#include "gaugefix.h"
#include "utilities.h"
#include "algebra.h"

//////////////////////////auxiliars
void setidLatticeGaugeSU2(LatticeGaugeSU2 * gauge_transform){
	int t,x,y,z;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					setidv(getg(gauge_transform,t,x,y,z));
				}
			}
		}
	}
}

void reunitLatticeGaugeSU2(LatticeGaugeSU2 * g){
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

double calcEps(LatticeLinkSU2 * lattice , LatticeGaugeSU2 * g){
	int t,x,y,z,mi,i;
	double sum = 0e0;
	double * auxU = malloc(sizeof(double)*4);
	int ami[4];
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			setUnitVector(&ami[0],mi);

			hermcv(auxU , getg(g, getStepT(t,+ami[0]) , getStep(x,+ami[1]) , getStep(y,+ami[2]) , getStep(z,+ami[3])) );
			mmprodv(auxU , getLink(lattice,t,x,y,z,mi) , auxU);
			mmprodv(auxU , getg(g,t,x,y,z) , auxU);
			sum += tracev(auxU);
		}
	}}}}
	sum /= (2.0*4*totalV);
	return 1-sum;
}

double applyFix(LatticeLinkSU2 * latticeOut, LatticeLinkSU2 * latticeIn , LatticeGaugeSU2 * g){
	clock_t stime = clock();
	double * auxU = malloc(sizeof(double)*4);
	int t,x,y,z,mi,i,ami[4];
	for(t=0 ; t<latticeIn->Nt ; t++){
	for(x=0 ; x<latticeIn->N ; x++){
	for(y=0 ; y<latticeIn->N ; y++){
	for(z=0 ; z<latticeIn->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			setUnitVector(&ami[0],mi);
			hermcv(auxU , getg(g, getStepT(t,+ami[0]) , getStep(x,+ami[1]) , getStep(y,+ami[2]) , getStep(z,+ami[3])) );
			mmprodv(auxU , getLink(latticeIn,t,x,y,z,mi) , auxU);
			mmprodv(getLink(latticeOut,t,x,y,z,mi) , getg(g,t,x,y,z) , auxU);
		}
	}}}}
	free(auxU);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

void copyg(LatticeGaugeSU2 * gout, LatticeGaugeSU2 * gin){
	int i;
	//check dimensions
	if((gout->N != gin->N) || (gout->Nt != gin->Nt))
		printf("\nWARNING! Copying vectors of different dimensions\n");
	for(i=0 ; i<gout->dimg ; i++){
		gout->g[i] = gin->g[i];
	}
}

double su2distance(LatticeLinkSU2 * lattice, LatticeLinkSU2 * W){
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
			cmprodv(aux1,-1,getLink(W,t,x,y,z,mi));
			sumv(aux1, getLink(lattice,t,x,y,z,mi) , aux1);
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

double los_alamos(LatticeLinkSU2 * lattice , LatticeGaugeSU2 * g , double * e2){
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
			setUnitVector(&ami[0],mi);

			hermcv( aux1 , getLink(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getLink(lattice,t,x,y,z,mi) , aux2);
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

double cornell(LatticeLinkSU2 * lattice , LatticeGaugeSU2 * g, double * e2 , double alpha){
	clock_t stime = clock();
	int t,x,y,z,mi,i;
	int ami[4];
	double auxN,auxT;
	double * id = malloc(sizeof(double)*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	*e2 = 0e0;
	setidv(id);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setUnitVector(&ami[0],mi);

			hermcv( aux1 , getLink(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getLink(lattice,t,x,y,z,mi) , aux2);
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

		*e2 += auxN*auxN*(1e0 - auxT*auxT/4e0);
	}}}}

	*e2 /= (double)(lattice->totalV);
	free(h);
	free(aux1);
	free(aux2);
	free(id);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

double overrelaxation(LatticeLinkSU2 * lattice , LatticeGaugeSU2 * g, double * e2 , double omega){
	clock_t stime = clock();
	int t,x,y,z,mi,i;
	int ami[4];
	double auxN,auxT;
	double * id = malloc(sizeof(double)*4);
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	*e2 = 0.E0;
	setidv(id);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setUnitVector(&ami[0],mi);

			hermcv( aux1 , getLink(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getLink(lattice,t,x,y,z,mi) , aux2);
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

		*e2 += auxN*auxN*(1e0 - auxT*auxT/4e0);
	}}}}
	*e2 /= (double)(lattice->totalV);

	free(id);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC );
}

double stochastic(LatticeLinkSU2 * lattice , LatticeGaugeSU2 * g , double * e2, double p){
	clock_t stime = clock();
	int t,x,y,z,mi,i,ami[4];
	double auxN,auxT;
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	*e2 = 0e0;
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setUnitVector(&ami[0],mi);

			hermcv( aux1 , getLink(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			mmprodv(aux1 , aux1 , aux2);
			sumv(h , h , aux1);

			hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			mmprodv(aux1 , getLink(lattice,t,x,y,z,mi) , aux2);
			sumv(h , h , aux1);
		}

		auxN = sqrt(detv(h));						//auxN <- sqrt(|h|) = N
		cmprodv(h , 1e0/auxN , h);					//h <- h/N = ~h
		mmprodv(aux1 , getg(g,t,x,y,z) , h);		//aux1 <- g . ~h = ~w
		auxT = tracev(aux1);						//aux1 <- tr(~w) = T

		if( ran0(global_seed) < p ){
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
			printf("\n inside stoch e2 = %.2e ; auxN = %.3e ; auxT = %.3e\n", *e2, auxN, auxT);
			getchar();
		}
	}}}}

	*e2 /= (double)(lattice->totalV);
	free(h);
	free(aux1);
	free(aux2);
	return( ((double)(clock() - stime))/CLOCKS_PER_SEC);
}

//AUTOMATIZATION
double fixLatticeStoch(LatticeLinkSU2 * lout, LatticeLinkSU2 * lin, LatticeGaugeSU2 * g, double p_stoch, double e2tol){
	clock_t stime = clock();
	printf("\nSO Landau Gauge Fixing initiated. Expected tolerance e2tol = %.0e\n",e2tol);
	//double * g 			= malloc(sizeof(double)*totalV*4);
	//setidLatticeGaugeSU2(g);
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

double fixLatticeLosalamos(LatticeLinkSU2 * lout, LatticeLinkSU2 * lin, LatticeGaugeSU2 * g , double e2tol){
	clock_t stime = clock();
	printf("\n Los Alamos gauge fixing initiated. Expected tolerance e2tol = %.0e\n",e2tol);
	//setidLatticeGaugeSU2(g);
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

double fixLatticeOver(LatticeLinkSU2 * lout, LatticeLinkSU2 * lin, LatticeGaugeSU2 * g, double omega, double e2tol){
	clock_t stime = clock();
	printf("\nSO Landau Gauge Fixing initiated. Expected tolerance e2tol = %.0e\n",e2tol);
	//double * g 			= malloc(sizeof(double)*totalV*4);
	//setidLatticeGaugeSU2(g);
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

double calce2(LatticeLinkSU2 * lattice){
	int t,x,y,z,mi,i;
	double sum = 0;
	double * h = malloc(sizeof(double)*4);
	double * aux1 = malloc(sizeof(double)*4);
	double * aux2 = malloc(sizeof(double)*4);
	double e2 = 0e0,auxN,auxT;
	int ami[4];
	LatticeGaugeSU2 * g = newLatticeGaugeSU2(lattice->N,lattice->Nt);
	setidLatticeGaugeSU2(g);
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		setzerov(h);
		for(mi=0 ; mi<4 ; mi++){
			setUnitVector(&ami[0],mi);

			hermcv(aux1, getLink(lattice,t,x,y,z,mi));
			sumv(aux1, aux1, getLink(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi));
			sumv(h,h,aux1);

			hermcv(aux2, getLink(lattice , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi));
			sumv(aux2,aux2,getLink(lattice,t,x,y,z,mi));
			cmprodv(aux2, -1, aux2);
			sumv(h,h,aux2);
			//hermcv( aux1 , getLink(U , getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3]) ,mi) );
			//hermcv( aux2 , getg(g, getStepT(t,-ami[0]) , getStep(x,-ami[1]) , getStep(y,-ami[2]) , getStep(z,-ami[3])) );
			//mmprodv(aux1 , aux1 , aux2);
			//sumv(h , h , aux1);

			//hermcv( aux2 , getg(g, getStepT(t,ami[0]) , getStep(x,ami[1]) , getStep(y,ami[2]) , getStep(z,ami[3]) ) );
			//mmprodv(aux1 , getLink(U,t,x,y,z,mi) , aux2);
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
	freeLatticeGaugeSU2(g);
	free(h);
	free(aux1);
	free(aux2);
	return e2/(totalV);
}

double calibrate_over(LatticeLinkSU2 * lattice , double tol){
	// the parameter omega from overrelaxation is supposed to obey 1<omega<2
	int i;
	double e2;
	LatticeGaugeSU2 * g = newLatticeGaugeSU2(lattice->N,lattice->Nt);
	double step_size = 0.1;
	double best_omega = 1.05E0;
	double omega_temp = 1.05E0;
	double best_path = 1e20;
	while(omega_temp < 2.0){
		setidLatticeGaugeSU2(g);
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
	freeLatticeGaugeSU2(g);
	return(best_omega);
}

double calibrate_cornell(LatticeLinkSU2 * lattice , double tol, double step_size, double lower_alpha, double upper_alpha){
	//the alpha parameter from the Cornell method is only supposed to be small and positive
	//so its calibration is more difficult
	int i,j;
	int n_steps = (int)((upper_alpha-lower_alpha)/step_size);
	double e2;
	double best_alpha;
	double temp_alpha;
	double best_path = 1e20;
	LatticeGaugeSU2 * g = newLatticeGaugeSU2(lattice->N,lattice->Nt);
	printf("\nCornell calibration initiated:\n");
	for(j=0 ; j < n_steps ; j++){
		setidLatticeGaugeSU2(g);
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
	freeLatticeGaugeSU2(g);
	return(best_alpha);
}

double calibrate_stoc(LatticeLinkSU2 * lattice , double tol){
	// the p parameter of stochastic overrelaxation is supposed to obey 0<p<1
	int i;
	double e2;
	double best_p = 0.85;
	double temp_p = best_p;
	double best_path = 1e20;
	double step_size = 0.075;
	LatticeGaugeSU2 * g = newLatticeGaugeSU2(lattice->N,lattice->Nt);
	LatticeLinkSU2 * latticetemp = newLatticeLinkSU2(lattice->N,lattice->Nt);
	printf("\nStochastic Overrelaxation calibration initiated:\n");
	while(temp_p > 0e0 ){
		setidLatticeGaugeSU2(g);
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
	freeLatticeGaugeSU2(g);
	freeLatticeLinkSU2(latticetemp);
	return(best_p);
}

//In Test
double calcFunctional(LatticeLinkSU2 * lattice){
	int t,x,y,z,mi,i;
	double sum = 0e0;
	int ami[4];
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			sum += tracev( getLink(lattice, t, x, y, z, mi) );
			//printf("\n tr u = %lf", tracev(getLink(U, t, x, y, z, mi)) );
		}
	}}}}
	sum /= (8.0*lattice->totalV);
	return 1-sum;
}

