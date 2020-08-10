#include "inverters.h"

#include "utilities.h"
#include "algebra.h"

const int ITERATION_TOL = 5000;

//CG and PM with real vectors
double rfpapply(LatticeColorVectorReal * Mvout , LatticeLinkSU2 * lattice , LatticeColorVectorReal * vin){
	//dim should be = N_colors * N^3 * Nt
	//the color index = a is fixed
	//ATTENTION: a is 1,2,3 or 0,1,2?
	clock_t dtime = clock();
	int a,b,c,t,x,y,z,mi;
	int ami[4];
	LatticeColorVectorReal * Mvouttemp = newLatticeColorVectorReal(vin->N,vin->Nt);
	// double * Mvouttemp = malloc(sizeof(double)*dim);

	setzeroLatticeColorVectorReal(Mvouttemp);

	for(t=0;t<lattice->Nt;t++){
	for(x=0;x<lattice->N;x++){
	for(y=0;y<lattice->N;y++){
	for(z=0;z<lattice->N;z++){
		for(a=0;a<3;a++){
			for(mi=0;mi<4;mi++){

				setzeroqvector(ami);
				ami[mi]=1;


				*getLatticeColorVectorReal(Mvouttemp,a,t,x,y,z) +=

				getLink(lattice,t,x,y,z,mi)[0]
				*(
				*getLatticeColorVectorReal(vin,a,t,x,y,z)
				- *getLatticeColorVectorReal(vin,a,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
				)

				+

				getLink(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[0]
				*(
				*getLatticeColorVectorReal(vin,a,t,x,y,z)
				-*getLatticeColorVectorReal(vin,a,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
				);


				//sum in b,c
				for(b=0;b<3;b++){
					for(c=0;c<3;c++){
						if((a!=b) && (b!=c) && (a!=c)){
							*getLatticeColorVectorReal(Mvouttemp,a,t,x,y,z) +=
							aseps(a,b,c)*
							(

								getLink(lattice,t,x,y,z,mi)[b+1]*
								(
									*getLatticeColorVectorReal(vin,c,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
									-
									0//*getelvr(vin,c,t,x,y,z)
								)

								-

								getLink(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[b+1]*
								(
									*getLatticeColorVectorReal(vin,c,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
									-
									0//*getelvr(vin,c,t,x,y,z)
								)

							);
						}
					}
				}
			}
		}
	}}}}
	copyLatticeColorVectorReal(Mvout,Mvouttemp);
	freeLatticeColorVectorReal(Mvouttemp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double rfpapplyOrthogonal(LatticeColorVectorReal * Mvout , LatticeLinkSU2 * lattice , LatticeColorVectorReal * vin){
	//dim should be = N_colors * N^4
	//the color index = a is fixed
	//ATTENTION: a is 1,2,3 or 0,1,2?
	clock_t dtime = clock();
	int a,b,c,t,x,y,z,mi;
	int ami[4];
	LatticeColorVectorReal * Mvouttemp = newLatticeColorVectorReal(vin->N,vin->Nt);

	setzeroLatticeColorVectorReal(Mvouttemp);

	for(t=0;t<Nt;t++){
	for(x=0;x<N;x++){
	for(y=0;y<N;y++){
	for(z=0;z<N;z++){
		for(a=0;a<3;a++){
			for(mi=0;mi<4;mi++){

				setzeroqvector(ami);
				ami[mi]=1;


				*getLatticeColorVectorReal(Mvouttemp,a,t,x,y,z) +=

				getLink(lattice,t,x,y,z,mi)[0]
				*(
				*getLatticeColorVectorReal(vin,a,t,x,y,z)
				-*getLatticeColorVectorReal(vin,a,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
				)

				+

				getLink(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[0]
				*(
				*getLatticeColorVectorReal(vin,a,t,x,y,z)
				-*getLatticeColorVectorReal(vin,a,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
				);


				//sum in b,c

				for(b=0;b<3;b++){
					for(c=0;c<3;c++){
						if((a!=b) && (b!=c) && (a!=c)){
							*getLatticeColorVectorReal(Mvouttemp,a,t,x,y,z) +=
							aseps(a,b,c)*
							(

								getLink(lattice,t,x,y,z,mi)[b+1]*
								(
									*getLatticeColorVectorReal(vin,c,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
									-
									*getLatticeColorVectorReal(vin,c,t,x,y,z)
								)

								-

								getLink(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[b+1]*
								(
									*getLatticeColorVectorReal(vin,c,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
									-
									*getLatticeColorVectorReal(vin,c,t,x,y,z)
								)

							);
						}
					}
				}
			}
		}
	}}}}
	copyLatticeColorVectorReal(Mvout,Mvouttemp);
	freeLatticeColorVectorReal(Mvouttemp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double rfpcgmr(LatticeColorVectorReal * x , LatticeLinkSU2 * lattice , LatticeColorVectorReal * b,  int dim, double tol, int * iCG, double * r_CG){
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a hermitian (dim x dim) matrix (real and symmetric in this case)
	//A is the FP matrix defined from the lattice links applied to the source b
	//in the usual context dim = N_colors * N^4 ==> the vector will have a super-index (a,x)
	//returns last r_CG value

	clock_t dtime = clock();
	int i;
	(*iCG)=0;
	double alpha, rnew,rold;;
	LatticeColorVectorReal * auxv = newLatticeColorVectorReal(lattice->N,lattice->Nt);
	LatticeColorVectorReal * r = newLatticeColorVectorReal(lattice->N,lattice->Nt);
	LatticeColorVectorReal * p = newLatticeColorVectorReal(lattice->N,lattice->Nt);

	//chute inicial passado em x
	//definições iniciais:

	rfpapply(auxv , lattice , x);		//mvprodr(auxv , A , x , dim);

	cprodLatticeColorVectorReal(auxv , -1e0 , auxv);
	sumLatticeColorVectorReal(r , b , auxv);				//r <- b - A.x

	copyLatticeColorVectorReal(p , r);						//p <- r

	rold = inprodLatticeColorVectorReal(r,r);

	// printf("\n Iter = %d ; |r_CG| = %.0e ",*iCG,sqrt(rnew));
	fflush(stdout);

	while(1){														//for i = 1:length(b)
		(*iCG)++;

		rfpapply(auxv,lattice,p);			//Ap = A * p;


		alpha = rold/inprodLatticeColorVectorReal(p,auxv);		//alpha = rsold / (p' * Ap);

		cprodLatticeColorVectorReal(auxv , -alpha , auxv);
		sumLatticeColorVectorReal(r , r , auxv);				//r = r - alpha * Ap;

		cprodLatticeColorVectorReal(auxv , alpha , p);
		sumLatticeColorVectorReal(x , x , auxv);				//x = x + alpha * p;

		rnew = inprodLatticeColorVectorReal(r,r);				//rsnew = r' * r;

		// printf("\r Iter = %d ; |r_CG| = %.0e ",(*iCG),sqrt(rnew));
		fflush(stdout);

		if( (sqrt(rnew) < tol) || (sqrt(rnew) > 1e8) || *iCG > ITERATION_TOL)					//if sqrt(rsnew) < tol
			break;								//break

		cprodLatticeColorVectorReal(auxv , rnew/rold , p);
		sumLatticeColorVectorReal(p , r , auxv);				//p = r + (rsnew / rsold) * p;

		rold = rnew;							//rsold = rsnew;
	}

	freeLatticeColorVectorReal(r);
	freeLatticeColorVectorReal(p);
	freeLatticeColorVectorReal(auxv);

	*r_CG = rnew;

	return((double)(clock()-dtime)/CLOCKS_PER_SEC);
}

void rgetOrthogonalSpace(LatticeColorVectorReal * c, LatticeColorVectorReal * b, LatticeLinkSU2 * lattice, double cg_tol, int * iCG, double * r_CG){
	LatticeColorVectorReal * a = newLatticeColorVectorReal(lattice->N,lattice->Nt);

	//getting the null space out of M
	rfpapply(a , lattice , b);
	//solve CG
	copyLatticeColorVectorReal(c,a);
	//double pmin[4] = {0e0,0e0,0e0,1e0/N};
	//rsetplanewave(c, pmin ,0);
	//setonevc(c,colorV);
	rfpcgmr(c, lattice , a, colorV, cg_tol, iCG, r_CG);

	freeLatticeColorVectorReal(a);
}

int rsmallest_eigen_cg(LatticeLinkSU2 * lattice, double * lambda1, LatticeColorVectorReal * eigen_out
	 	, double eigen_tol , LatticeColorVectorReal * vguess, double * r_cg, int * iCG, double cg_tol){
	
	clock_t dtime = clock();
	double norm;
	double eigen_v = 0e0;
	double eigen_vp = 0e0;

	LatticeColorVectorReal * X = newLatticeColorVectorReal(lattice->N,lattice->Nt);
	LatticeColorVectorReal * Y = newLatticeColorVectorReal(lattice->N,lattice->Nt);
	LatticeColorVectorReal * Xp = newLatticeColorVectorReal(lattice->N,lattice->Nt);

	*r_cg = 0e0;
	int i_PM = 0;
	//double pmin[4] = {0e0,0e0,0e0,1e0/N};
	//rsetplanewave(X, pmin ,2);
	//rsetplanewaveallcolors(X, pmin);
	//setonevr(X,colorV);
	//setrandomvr(X, colorV, _seed);

	copyLatticeColorVectorReal(X,vguess);											//vguess is some guess to the eigenvector with the smallest eigenvalue

	//Power method

	//rgetOrthogonalSpace(X,X,lattice);
	printf("\nPM: Iter=%d | lambda_1=%e |  lambda1_norm=%e",i_PM,1e0/eigen_vp,norm);
	fflush(stdout);

	do{
		//if(i%3 == 0)
		//printf(" \n->-> Average projection on plane-wave = %e\n", compareVectorToPW(X, 0));
		//if(compareVectorToPW(X, 0) > 1e-10)
		rgetOrthogonalSpace(X,X,lattice,cg_tol, iCG, r_cg);
		if(*r_cg > 1e2 || *iCG>ITERATION_TOL){
			printf("\n-----> r_CG DIVERGED in PM (orthogonalization)!\n");
			fflush(stdout);
			//getchar();
			break;
		}
		//getchar();
		//rfpapply(X , lattice , X , colorV);
		copyLatticeColorVectorReal(Xp,X);

		rfpcgmr(Xp , lattice , X , colorV, cg_tol, iCG, r_cg);
		if(*r_cg > 1e2 || *iCG>ITERATION_TOL){
			printf("\n-----> r_CG DIVERGED in PM\n");
			fflush(stdout);
			//getchar();
			break;
		}

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodLatticeColorVectorReal(X,Xp)/inprodLatticeColorVectorReal(X,X);
		copyLatticeColorVectorReal(X,Xp);
		reunitLatticeColorVectorReal(X);

		i_PM++;
		// eigen_vp = (eigen_vp);
		rfpapply(Xp,lattice,X);
		cprodLatticeColorVectorReal(Y,-1e0/eigen_vp,X);
		sumLatticeColorVectorReal(Xp,Xp,Y);
		norm = normLatticeColorVectorReal(Xp);
		printf("\nPM: Iter=%d | lambda_1=%e |  lambda1_norm=%e",i_PM,1e0/eigen_vp,norm);
		fflush(stdout);

		// cprodvr(X, getSignal(min(X,colorV)), X, colorV);		//Attilio's trick
	}while( norm > eigen_tol);	//this way the eigenvalue equation is directly verified, but is slower

	//	printf("\n ev = %lf | evp = %lf | abs() = %lf | tol = %lf \n",1.0/eigen_v,1.0/eigen_vp,fabs(1.0/eigen_v - 1.0/eigen_vp),eigen_tol);
	//}while( fabs(1.0/eigen_v - 1.0/eigen_vp) > 1e-7 ); 	//this way is faster but the tolerance must be smaller

	copyLatticeColorVectorReal(eigen_out, X);

	freeLatticeColorVectorReal(X);
	freeLatticeColorVectorReal(Xp);
	freeLatticeColorVectorReal(Y);
	*lambda1 = 1e0/eigen_vp;
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

//CG and PM with complex vector ; still to wrap in structs
double fpapply(LatticeColorVectorComplex * Mvout , LatticeLinkSU2 * lattice , LatticeColorVectorComplex * vin){
	clock_t dtime = clock();
	//dim should be = N_colors * N^4
	//the color index = a is fixed
	//ATTENTION: is a=1,2,3 or 0,1,2?
	int a,b,c,t,x,y,z,mi;
	int ami[4];
	LatticeColorVectorComplex * Mvouttemp = newLatticeColorVectorComplex(lattice->N,lattice->Nt);

	setzeroLatticeColorVectorComplex(Mvouttemp);

	for(a=0;a<3;a++){
		for(t=0;t<lattice->Nt;t++){
		for(x=0;x<lattice->N;x++){
		for(y=0;y<lattice->N;y++){
		for(z=0;z<lattice->N;z++){
			for(mi=0;mi<4;mi++){
				setUnitVector(&ami[0],mi);

				*getLatticeColorVectorComplex(Mvouttemp,a,t,x,y,z) +=

				getLink(lattice,t,x,y,z,mi)[0]*
				(
				*getLatticeColorVectorComplex(vin,a,t,x,y,z)
				-
				*getLatticeColorVectorComplex(vin,a,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
				)
				+
				getLink(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[0]*
				(
				*getLatticeColorVectorComplex(vin,a,t,x,y,z)
				-
				*getLatticeColorVectorComplex(vin,a,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
				);

				//sum in b,c
				for(b=0;b<3;b++){
					for(c=0;c<3;c++){
						if((a!=b) && (b!=c)){
							*getLatticeColorVectorComplex(Mvouttemp,a,t,x,y,z) +=
							aseps(a,b,c)*
							(

							getLink(lattice,t,x,y,z,mi)[b+1]*
							(
								*getLatticeColorVectorComplex(vin,c,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
								-
								0//*getLatticeColorVectorComplex(vin,c,t,x,y,z)
							)

							-

							getLink(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[b+1]*
							(
								*getLatticeColorVectorComplex(vin,c,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
								-
								0//*getLatticeColorVectorComplex(vin,c,t,x,y,z)
							)
							);
						}
					}
				}
			}
		}}}}
	}

	copyLatticeColorVectorComplex(Mvout,Mvouttemp);
	freeLatticeColorVectorComplex(Mvouttemp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double fpcgmr(LatticeColorVectorComplex * x , LatticeLinkSU2 * lattice , LatticeColorVectorComplex * b){
	clock_t dtime = clock();
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a hermitian (dim x dim) matrix (real and symmetric in this case)
	//A is the FP matrix defined from the lattice links applied to the source b
	//in the usual context dim = N_colors * N^4 ==> the vector will have a super-index (a,x)

	double tol = 1e-10;
	int i,icont = 0;
	double complex alpha, rnew,rold;

	LatticeColorVectorComplex * auxv = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * r = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * p = newLatticeColorVectorComplex(lattice->N,lattice->Nt);

	//chute inicial passado em x
	//definições iniciais:

	fpapply(auxv , lattice , x);		//mvprodr(auxv , A , x , dim);

	//fpapply(auxv , lattice , auxv  ,dim);		//mvprodr(auxv , A , x , dim);
	//fpapply(b , lattice , b  ,dim);		//mvprodr(auxv , A , x , dim);

	cprodLatticeColorVectorComplex(auxv , -1e0 , auxv);
	sumLatticeColorVectorComplex(r , b , auxv );				//r <- b - A.x

	copyLatticeColorVectorComplex(p , r);						//p <- r

	rold = inprodLatticeColorVectorComplex(r,r);

	// printf("\n|r_CG| = %.0e ",sqrt(rnew));
	fflush(stdout);

	while(1){														//for i = 1:length(b)
		icont++;

		fpapply(auxv,lattice,p);			//Ap = A * p;

		alpha = rold/inprodLatticeColorVectorComplex(p,auxv);		//alpha = rsold / (p' * Ap);

		cprodLatticeColorVectorComplex(auxv , -alpha , auxv);
		sumLatticeColorVectorComplex(r , r , auxv);				//r = r - alpha * Ap;

		cprodLatticeColorVectorComplex(auxv , alpha , p);
		sumLatticeColorVectorComplex(x , x , auxv);				//x = x + alpha * p;

		rnew = inprodLatticeColorVectorComplex(r,r);				//rsnew = r' * r;

		// printf("\r|r_CG| = %.0e ",sqrt(rnew));
		fflush(stdout);

		if(sqrt(rnew) < tol)					//if sqrt(rsnew) < tol
			break;								//break

		cprodLatticeColorVectorComplex(auxv , rnew/rold , p);
		sumLatticeColorVectorComplex(p , r , auxv);				//p = r + (rsnew / rsold) * p;

		rold = rnew;							//rsold = rsnew;
	}

	freeLatticeColorVectorComplex(r);
	freeLatticeColorVectorComplex(p);
	freeLatticeColorVectorComplex(auxv);

	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

void getOrthogonalSpace(LatticeColorVectorComplex * c, LatticeColorVectorComplex * b, LatticeLinkSU2 * lattice){
	LatticeColorVectorComplex * a = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	//getting the null space out of M
	fpapply(a , lattice , b);
	//solve CG
	copyLatticeColorVectorComplex(c,a);
	//setonevc(c,colorV);
	fpcgmr(c, lattice , a);

	freeLatticeColorVectorComplex(a);
}

double complex smallest_eigen_cg(LatticeLinkSU2 * lattice, double * lambda1, LatticeColorVectorComplex * eigen_out){
	clock_t dtime = clock();
	int i;
	double norm;
	double eigen_tol = 1e-5;
	double complex eigen_v = 0e0 + I*0e0;
	double complex eigen_vp = 0e0 + I*0e0;

	LatticeColorVectorComplex * X = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * Y = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * Xp = newLatticeColorVectorComplex(lattice->N,lattice->Nt);

	double pmin[4] = {0e0,0e0,0e0,1e0/N};
	setplanewaveallcolorsLatticeColorVectorComplex(X, pmin);

	//Power method
	i=0;
	printf("\nIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);
	do{
		getOrthogonalSpace(X,X,lattice);
		copyLatticeColorVectorComplex(Xp,X);
		fpcgmr(Xp , lattice , X);

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodLatticeColorVectorComplex(X,Xp)/inprodLatticeColorVectorComplex(X,X);
		copyLatticeColorVectorComplex(X,Xp);
		reunitLatticeColorVectorComplex(X);

		i++;

		fpapply(Xp,lattice,X);
		cprodLatticeColorVectorComplex(Y,-1e0/eigen_vp,X);
		sumLatticeColorVectorComplex(Xp, Xp,Y);

		norm = normLatticeColorVectorComplex(Xp);

		printf("\rIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);

	}while( norm > eigen_tol );
	copyLatticeColorVectorComplex(eigen_out, X);

	freeLatticeColorVectorComplex(X);
	freeLatticeColorVectorComplex(Xp);
	freeLatticeColorVectorComplex(Y);
	*lambda1 = 1e0/creal(eigen_vp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

//biCGstab with complex vectors
void fpbicgstabmr(LatticeColorVectorComplex * x , LatticeLinkSU2 * lattice , LatticeColorVectorComplex * b, double cg_tol, double * r_cg, int * icg){
	clock_t dtime = clock();
	//algorithm based along the lines of Gattringer's book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a complex (dim x dim) matrix
	//A is the FP matrix defined from the lattice links
	//in the usual context dim = N_colors * N^4 ==> the vector will have a super-index (a,x)
	double tol = cg_tol;
	int i = 0;
	*icg=0;
	double complex alpha, rho1, rho2 , omega , beta;

	LatticeColorVectorComplex * aux1 = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * aux2 = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * rtilde = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * r = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * v = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * s = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * t = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * p = newLatticeColorVectorComplex(lattice->N,lattice->Nt);

	//chute inicial passado em x
	//definições iniciais:

	fpapply(aux1 , lattice , x);		//aux1 <- M.x

	cprodLatticeColorVectorComplex(aux1 , -1e0 , aux1);
	sumLatticeColorVectorComplex(r , b , aux1);				//r <- b - A.x

	copyLatticeColorVectorComplex(rtilde, r);						//rtilde <- r

	printf("\n iter=%d |s_biCGstab| = %e", *icg, *r_cg);
	fflush(stdout);
	while(1){//for(i=0;i<dim;i++){				//for i = 1:length(b)
		(*icg)++;

		rho1 = inprodLatticeColorVectorComplex(rtilde, r);

		if( cabs(rho1) == 0){
			printf("\nrho=0 : BiCGStab failed.");
			break;
		}

		if(*icg==1){
			copyLatticeColorVectorComplex(p,r);
		}
		else{
			beta = alpha*rho1/(omega*rho2);

			cprodLatticeColorVectorComplex(aux1,-omega,v);
			sumLatticeColorVectorComplex(aux1,aux1,p);
			cprodLatticeColorVectorComplex(aux1,beta,aux1);
			sumLatticeColorVectorComplex(p,r,aux1);			//p <- r + beta.(p-omega.v)
		}

		fpapply(v,lattice,p);			//v <- A.p

		alpha = rho1/inprodLatticeColorVectorComplex(rtilde,v);	//alpha <- rho1/(rtilde,v)

		cprodLatticeColorVectorComplex(aux1,-alpha,v);
		sumLatticeColorVectorComplex(s,r,aux1);

		*r_cg = normLatticeColorVectorComplex(s);

		printf("\r iter=%d |s_biCGstab| = %e", *icg, *r_cg);
		fflush(stdout);

		if( *r_cg < tol || *r_cg > 1e2 || *icg>ITERATION_TOL){
			cprodLatticeColorVectorComplex(aux1,alpha,p);
			sumLatticeColorVectorComplex(x,x,aux1);			//x <- x + alpha.p
			break;
		}

		fpapply(t,lattice,s);			//t <- A.s

		omega = inprodLatticeColorVectorComplex(t,s)/inprodLatticeColorVectorComplex(t,t);	//omega <- (t,s)/(t,t)

		cprodLatticeColorVectorComplex(aux1,-omega,t);
		sumLatticeColorVectorComplex(r,s,aux1);				//r <- s - omega.t

		cprodLatticeColorVectorComplex(aux1,alpha,p);
		cprodLatticeColorVectorComplex(aux2,omega,s);
		sumLatticeColorVectorComplex(aux1,aux1,aux2);
		sumLatticeColorVectorComplex(x,x,aux1);				//x <- x + alpha.p + oemga.s

		rho2 = rho1;
	}

	freeLatticeColorVectorComplex(r);
	freeLatticeColorVectorComplex(rtilde);
	freeLatticeColorVectorComplex(v);
	freeLatticeColorVectorComplex(s);
	freeLatticeColorVectorComplex(t);
	freeLatticeColorVectorComplex(p);
	freeLatticeColorVectorComplex(aux1);
	freeLatticeColorVectorComplex(aux2);

	//return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

void bigetOrthogonalSpace(LatticeColorVectorComplex * c, LatticeColorVectorComplex * b, LatticeLinkSU2 * lattice, double cg_tol, double * r_cg, int * icg){
	LatticeColorVectorComplex * a = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	//getting the null space out of M
	fpapply(a , lattice , b);
	//solve CG
	copyLatticeColorVectorComplex(c,a);
	//setonevc(c,colorV);
	fpbicgstabmr(c, lattice , a, cg_tol, r_cg, icg);

	freeLatticeColorVectorComplex(a);
}

//R functional computation

double  FsecondDerivative(LatticeLinkSU2 * lattice, LatticeColorVectorComplex * eigen_vector_out){
	double lambda;
	smallest_eigen_cg(lattice,&lambda,eigen_vector_out);
	return lambda/(4e0*lattice->totalV);
}

double rFThirdFourthDerivative(LatticeLinkSU2 * lattice, double * ev_in, double * third, double * fourth){
	clock_t dtime = clock();
	int t,y,x,z,a,b,c,d,e;
	double aux1,aux2;
	*fourth = 0e0;
	*third = 0e0;
	int emi[4], mi;
	double gamu, ga, gamma2mu, gamma2, gammaDiv, gammaAmu, gammaMgamma;
	double ua, uamu;
	//double * Mgamma = malloc(sizeof(double)*colorV);

	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					for(mi=0;mi<4;mi++){
						setUnitVector(emi,mi);
						for(a=0;a<3;a++){
							for(b=0;b<3;b++){

								*third -=
								0.75*( pow(*getelvr(ev_in,a,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3])),2)
									    - pow(*getelvr(ev_in,a,t,x,y,z),2) )
									*( *getelvr(ev_in,b,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3]))
									    + *getelvr(ev_in,b,t,x,y,z) )
									*( getLink(lattice,t,x,y,z,mi)[b+1] );

								*third -=
								 (pow(*getelvr(ev_in,a,t,x,y,z),2) )
									*( *getelvr(ev_in,b,t,x,y,z) )
									*( getLink(lattice,t,x,y,z,mi)[b+1]
										- getLink(lattice,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]),mi)[b+1] );


								*fourth += 0.75*(
									pow(*getelvr(ev_in,a,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3])),2)
									-
									pow(*getelvr(ev_in,a,t,x,y,z),2)
								)
								*(
									pow(*getelvr(ev_in,b,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3])),2)
									-
									pow(*getelvr(ev_in,b,t,x,y,z),2)
								)
								*getLink(lattice,t,x,y,z,mi)[0];

								*fourth -= pow(*getelvr(ev_in,a,t,x,y,z),2)
											*( *getelvr(ev_in,b,t,x,y,z) )*
									(
										getLink(lattice,t,x,y,z,mi)[0]*(
											*getelvr(ev_in,b,t,x,y,z) -  *getelvr(ev_in,b,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3]))
										)

										+

										getLink(lattice,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]),mi)[0]*(
											*getelvr(ev_in,b,t,x,y,z) - *getelvr(ev_in,b,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]))
										)
									);

								for(d=0;d<3;d++){
									for(e=0;e<3;e++){
										if((b!=d) && (d!=e) && (b!=e)){
											*fourth -= pow(*getelvr(ev_in,a,t,x,y,z),2)*( *getelvr(ev_in,b,t,x,y,z) )*
												(
													aseps(b,d,e)*(
														getLink(lattice,t,x,y,z,mi)[d+1]*( *getelvr(ev_in,e,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3])) )
														-
														getLink(lattice,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]),mi)[d+1]*(
															*getelvr(ev_in,e,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3])) )
													)
												);
										}
									}
								}

							}
						}
					}
				}
			}
		}
	}
	*third /= totalV;
	*fourth /= totalV;

	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

//OTHER STUFF
double complex largest_eigen_cg(LatticeLinkSU2 * lattice, double * lambda1, LatticeColorVectorComplex * eigen_out){
	clock_t dtime = clock();
	int i;
	double eigen_tol = 1e-8;
	double complex eigen_v = 0e0 + I*0e0;
	double complex eigen_vp = 0e0 + I*0e0;

	LatticeColorVectorComplex * X = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * Y = newLatticeColorVectorComplex(lattice->N,lattice->Nt);
	LatticeColorVectorComplex * Xp = newLatticeColorVectorComplex(lattice->N,lattice->Nt);

	double pmin[4] = {0e0,0e0,0e0,1e0/N};
	setplanewaveallcolorsLatticeColorVectorComplex(X, pmin);
	//setonevc(X,colorV);

	double norm;

	//Power method
	i=0;
	//printf("\nIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);
	do{

		//CG
		//setplanewaveallcolors(Xp, pmin);
		//setonevc(Xp,colorV);
		//setzerovc(Xp,colorV);
		fpapply(Xp , lattice , X);

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodLatticeColorVectorComplex(X,Xp)/inprodLatticeColorVectorComplex(X,X);
		copyLatticeColorVectorComplex(X,Xp);
		reunitLatticeColorVectorComplex(X);

		i++;

		fpapply(Xp,lattice,X);
		cprodLatticeColorVectorComplex(Y,-eigen_vp,X);
		sumLatticeColorVectorComplex(Xp, Xp,Y);

		norm = normLatticeColorVectorComplex(Xp);

		//printf("\rIter=%d | lambda_1=%e |  lambda1_norm=%e",i,creal(eigen_vp),norm);

	}while( norm > eigen_tol );
	copyLatticeColorVectorComplex(eigen_out, X);

	freeLatticeColorVectorComplex(X);
	freeLatticeColorVectorComplex(Xp);
	freeLatticeColorVectorComplex(Y);
	*lambda1 = 1e0/creal(eigen_vp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double verifyOrthogonalization(LatticeLinkSU2 * lattice, LatticeColorVectorReal * vin){
	LatticeColorVectorReal * aux = newLatticeColorVectorReal(lattice->N,lattice->Nt);
	LatticeColorVectorReal * constant_vector	= newLatticeColorVectorReal(lattice->N,lattice->Nt);

	rfpapply(aux, lattice, vin);

	setoneLatticeColorVectorReal(constant_vector);

	double result = inprodLatticeColorVectorReal(constant_vector,aux);

	freeLatticeColorVectorReal(constant_vector);
	freeLatticeColorVectorReal(aux);

	return(result);
}

