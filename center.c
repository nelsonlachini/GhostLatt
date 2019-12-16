#include "center.h"
#include "global.h"

#include "gaugefix.h"

double largest_eigen(double * matrix, int dim,  double * lambda1, double  * eigen_out){
	int i;
	double eigen_tol = 1e-10;
	double eigen_v = 0e0 , eigen_vp = 0e0;
	double * X = malloc(sizeof(double)*dim);
	double * Y = malloc(sizeof(double)*dim);
	double * Xp = malloc(sizeof(double)*dim);
	setonevr(X,dim);
	reunitvr(X,dim);
	clock_t dtime = clock();
	double norm;

	//printmr(matrix,dim);

	//Power method
	i=0;
	//printf("\nIter=%d | lambda_1=%e |  lambda1_norm=%e",i,eigen_vp,norm);
	do{
		//CG
		mvrprod(Xp , matrix , X , dim);

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodvr(X,Xp,dim)/inprodvr(X,X,dim);
		copyvr(X,Xp,dim);
		reunitvr(X,dim);

		i++;
		//printf("\rIter=%d | lambda_1=%e |  lambda1_norm=%e",i,eigen_vp,norm);

	}while( fabs(eigen_vp - eigen_v) > eigen_tol );
	copyvr(eigen_out, X, dim);

	free(X);
	free(Xp);
	free(Y);
	*lambda1 = eigen_vp;
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

void constructD(double * lattice, double * D, int t, int x, int y, int z){
	int mu, i, j, emu[4] = {0,0,0,0};
	setzeromr(D,4);

	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			/*
			if(i==0 && j==0){
				for(mu=0;mu<4;mu++){
					setquadv(&emu[0],mu);
					//D00 case
					*(D+4*i+j) += getU(lattice,t,x,y,z,mu)[0]*getU(lattice,t,x,y,z,mu)[0];
					*(D+4*i+j) += getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)[0]
					*getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)[0];
				}
			}
			else */
			if( ((i!=0) && (j!=0)) || ((i==0) && (j==0)) ){
				for(mu=0;mu<4;mu++){
					setquadv(&emu[0],mu);
					//Dij, i,j=1,2,3 cases or both 0
					*(D+4*i+j) += getU(lattice,t,x,y,z,mu)[i]*getU(lattice,t,x,y,z,mu)[j];
					*(D+4*i+j) += getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)[i]
					*getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)[j];
				}
			}
			else{
				for(mu=0;mu<4;mu++){
					setquadv(&emu[0],mu);
					//D0i or Di0 cases
					*(D+4*i+j) += -getU(lattice,t,x,y,z,mu)[0]*getU(lattice,t,x,y,z,mu)[i+j];	//attention for the change of sign due to notation!
					*(D+4*i+j) += getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)[0]
					*getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)[i+j];
				}
			}
		}
	}
	//printf("\nD:");
	//printmr(D,4);
	//getchar();
}

double directMCGOverrelaxationSweep(double * lattice, double * g, double omega){
	int mu, t, x, y, z, emu[4] = {0,0,0,0};
	double D[4][4], lambda_max;
	double * aux1 = malloc(sizeof(double)*4);
	double deltaR=0;
	double norm, phi;
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					constructD(lattice, &D[0][0], t, x, y, z);
					largest_eigen(&D[0][0], 4,  &lambda_max,  getg(g,t,x,y,z));

					//overrelaxtion g^omega = Proj(unit + omega*(g-unit))
					cmprodv(getg(g,t,x,y,z), omega, getg(g,t,x,y,z));
					getg(g,t,x,y,z)[0] += 1 - omega;
					reunitv(getg(g,t,x,y,z));

					for(mu=0;mu<4;mu++){
						setquadv(&emu[0],mu);

						mmprodv(getU(lattice,t,x,y,z,mu) , getg(g,t,x,y,z) , getU(lattice,t,x,y,z,mu));

						hermcv(aux1, getg(g,t,x,y,z));
						mmprodv(getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)
						, getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu) , aux1);
					}

				}
			}
		}
	}
	free(aux1);
	return deltaR;
}

double directMCGStochasticSweep(double * lattice, double * g, double p_stoch){
	int mu, t, x, y, z, emu[4] = {0,0,0,0};
	double D[4][4], lambda_max;
	double * aux1 = malloc(sizeof(double)*4);
	double deltaR=0;
	double norm, phi;
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					constructD(lattice, &D[0][0], t, x, y, z);
					largest_eigen(&D[0][0], 4,  &lambda_max,  getg(g,t,x,y,z));


					if(ran0(global_seed) < p_stoch){
						//overrelaxtion g^omega = Proj(unit + omega*(g-unit))
						cmprodv(getg(g,t,x,y,z), 1.7, getg(g,t,x,y,z));
						getg(g,t,x,y,z)[0] += 1 - 1.7;
						reunitv(getg(g,t,x,y,z));
					}


					for(mu=0;mu<4;mu++){
						setquadv(&emu[0],mu);

						mmprodv(getU(lattice,t,x,y,z,mu) , getg(g,t,x,y,z) , getU(lattice,t,x,y,z,mu));

						hermcv(aux1, getg(g,t,x,y,z));
						mmprodv(getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu)
						, getU(lattice,getStepT(t,-emu[0]),getStep(x,-emu[1]),getStep(y,-emu[2]),getStep(z,-emu[3]),mu) , aux1);
					}

				}
			}
		}
	}
	free(aux1);
	return deltaR;
}

double calcR(double * lattice){
	double R=0;
	int t, x, y, z, mu;
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					for(mu=0;mu<4;mu++){
						R += getU(lattice,t,x,y,z,mu)[0]*getU(lattice,t,x,y,z,mu)[0];
					}
				}
			}
		}
	}
	return R;
}

void centerProjection(double * vortexlattice, double * MCGlattice){
	//center projection
	setzerovr(vortexlattice,dimLattice);
	int contv=0;
	int contnv=0;
	int t, x, y, z, mu;
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					for(mu=0;mu<4;mu++){
						getU(vortexlattice,t,x,y,z,mu)[0] = signal(getU(MCGlattice,t,x,y,z,mu)[0]);
						if(getU(vortexlattice,t,x,y,z,mu)[0] == 1){
							contnv++;
						}
						else if(getU(vortexlattice,t,x,y,z,mu)[0] == -1){
							contv++;
						}
					}
				}
			}
		}
	}
	//printf("\n Vortex count: %d vs %d",contnv,contv);
}

void centerRemoval(double * MCGlattice, double * vortexlattice){
	//vortex removal
	int t, x, y, z, mu;
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					for(mu=0;mu<4;mu++){
						//if( getU(vortexlattice,t,x,y,z,mu)[0] == -1 )
							//printm(getU(vortexlattice,t,x,y,z,mu));
							//printm(getU(MCGlattice,t,x,y,z,mu));
							mmprodv(getU(MCGlattice,t,x,y,z,mu), getU(vortexlattice,t,x,y,z,mu), getU(MCGlattice,t,x,y,z,mu));
							//printm(getU(MCGlattice,t,x,y,z,mu));
							//getchar();
							//cmprodv(getU(MCGlattice,t,x,y,z,mu),  getU(vortexlattice,t,x,y,z,mu)[0], getU(MCGlattice,t,x,y,z,mu));
						//printf("\n%d",signal(getU(fulllattice,t,x,y,z,mu)[0]));
					}
				}
			}
		}
	}
}

double MCGOver(double * outlattice, double * inlattice, /* double * g, */ double omega, double Rtol){
	clock_t dtime = clock();
	int i;
	double deltaR;
	double * g = malloc(sizeof(double)*totalV*4);
	initg(g);
	copyl(outlattice,inlattice);
	printf("\n MCG gauge-fixing:");
	printf("\n|deltaR|=%.2e ",fabs(deltaR));
	do{
		//printf("\n%lf",calcR(lattice));
		deltaR = calcR(outlattice);
		directMCGOverrelaxationSweep(outlattice, g, omega);
		deltaR -= calcR(outlattice);
		printf("\r|deltaR|=%.2e ",fabs(deltaR));
		fflush(stdout);
	}while( fabs(deltaR) > Rtol);
	//getchar();
	free(g);
	return((double)(clock()-dtime)/(CLOCKS_PER_SEC));
}

double MCGStoch(double * outlattice, double * inlattice, /* double * g, */ double p_stoch, double Rtol){
	clock_t dtime = clock();
	int i;
	double deltaR;
	double * g = malloc(sizeof(double)*totalV*4);
	initg(g);
	copyl(outlattice,inlattice);
	printf("\n MCG gauge-fixing:");
	printf("\n|deltaR|=%.2e ",fabs(deltaR));
	do{
		//printf("\n%lf",calcR(lattice));
		deltaR = calcR(outlattice);
		directMCGStochasticSweep(outlattice, g, p_stoch);
		deltaR -= calcR(outlattice);
		printf("\r|deltaR|=%.2e ",fabs(deltaR));
	}while( fabs(deltaR) > Rtol);
	//getchar();
	free(g);
	return((double)(clock()-dtime)/(CLOCKS_PER_SEC));
}

void centerDecomposeLattice(double * vortexlattice, double * vortexremovedlattice, double * fulllattice, double center_omega, double Rtol){
	MCGOver(vortexremovedlattice, fulllattice, center_omega, Rtol);
	centerProjection(vortexlattice, vortexremovedlattice);
	centerRemoval(vortexremovedlattice, vortexlattice);
}
