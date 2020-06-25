#include "inverters.h"

#include "utilities.h"
#include "algebra.h"

const int ITERATION_TOL = 5000;

//CG and PM with real vectors
double rfpapply(double * Mvout , LatticeSU2 * lattice , double * vin , int dim){
	//dim should be = N_colors * N^4
	//the color index = a is fixed
	//ATTENTION: a is 1,2,3 or 0,1,2?
	clock_t dtime = clock();
	int a,b,c,t,x,y,z,mi;
	int ami[4];
	double * Mvouttemp = malloc(sizeof(double)*dim);

	setzerovr(Mvouttemp,dim);

	for(t=0;t<lattice->Nt;t++){
	for(x=0;x<lattice->N;x++){
	for(y=0;y<lattice->N;y++){
	for(z=0;z<lattice->N;z++){
		for(a=0;a<3;a++){
			for(mi=0;mi<4;mi++){

				setzeroqvector(ami);
				ami[mi]=1;


				*getelvr(Mvouttemp,a,t,x,y,z) +=

				getU(lattice,t,x,y,z,mi)[0]
				*(
				*getelvr(vin,a,t,x,y,z)
				-*getelvr(vin,a,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
				)

				+

				getU(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[0]
				*(
				*getelvr(vin,a,t,x,y,z)
				-*getelvr(vin,a,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
				);


				//sum in b,c
				for(b=0;b<3;b++){
					for(c=0;c<3;c++){
						if((a!=b) && (b!=c) && (a!=c)){
							*getelvr(Mvouttemp,a,t,x,y,z) +=
							aseps(a,b,c)*
							(

								getU(lattice,t,x,y,z,mi)[b+1]*
								(
									*getelvr(vin,c,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
									-
									0//*getelvr(vin,c,t,x,y,z)
								)

								-

								getU(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[b+1]*
								(
									*getelvr(vin,c,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
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
	copyvr(Mvout,Mvouttemp,dim);
	free(Mvouttemp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double rfpapplyOrthogonal(double * Mvout , LatticeSU2 * lattice , double * vin , int dim){
	//dim should be = N_colors * N^4
	//the color index = a is fixed
	//ATTENTION: a is 1,2,3 or 0,1,2?
	clock_t dtime = clock();
	int a,b,c,t,x,y,z,mi;
	int ami[4];
	double * Mvouttemp = malloc(sizeof(double)*dim);

	setzerovr(Mvouttemp,dim);

	for(t=0;t<Nt;t++){
	for(x=0;x<N;x++){
	for(y=0;y<N;y++){
	for(z=0;z<N;z++){
		for(a=0;a<3;a++){
			for(mi=0;mi<4;mi++){

				setzeroqvector(ami);
				ami[mi]=1;


				*getelvr(Mvouttemp,a,t,x,y,z) +=

				getU(lattice,t,x,y,z,mi)[0]
				*(
				*getelvr(vin,a,t,x,y,z)
				-*getelvr(vin,a,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
				)

				+

				getU(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[0]
				*(
				*getelvr(vin,a,t,x,y,z)
				-*getelvr(vin,a,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
				);


				//sum in b,c

				for(b=0;b<3;b++){
					for(c=0;c<3;c++){
						if((a!=b) && (b!=c) && (a!=c)){
							*getelvr(Mvouttemp,a,t,x,y,z) +=
							aseps(a,b,c)*
							(

								getU(lattice,t,x,y,z,mi)[b+1]*
								(
									*getelvr(vin,c,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
									-
									*getelvr(vin,c,t,x,y,z)
								)

								-

								getU(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi)[b+1]*
								(
									*getelvr(vin,c,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
									-
									*getelvr(vin,c,t,x,y,z)
								)

							);
						}
					}
				}
			}
		}
	}}}}
	copyvr(Mvout,Mvouttemp,dim);
	free(Mvouttemp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double rfpcgmr(double * x , LatticeSU2 * lattice , double * b,  int dim, double tol, int * iCG, double * r_CG){
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
	double * auxv = malloc(sizeof(double)*dim);
	double * r = malloc(sizeof(double)*dim);
	double * p = malloc(sizeof(double)*dim);

	//chute inicial passado em x
	//definições iniciais:

	rfpapply(auxv , lattice , x  ,dim);		//mvprodr(auxv , A , x , dim);


	cprodvr(auxv , -1e0 , auxv , dim);
	sumvr(r , b , auxv , dim);				//r <- b - A.x

	copyvr(p , r,dim);						//p <- r

	rold = inprodvr(r,r,dim);

	printf("\n Iter = %d ; |r_CG| = %.0e ",*iCG,sqrt(rnew));
	fflush(stdout);

	while(1){														//for i = 1:length(b)
		(*iCG)++;

		rfpapply(auxv,lattice,p,dim);			//Ap = A * p;


		alpha = rold/inprodvr(p,auxv,dim);		//alpha = rsold / (p' * Ap);

		cprodvr(auxv , -alpha , auxv, dim);
		sumvr(r , r , auxv , dim);				//r = r - alpha * Ap;

		cprodvr(auxv , alpha , p ,dim);
		sumvr(x , x , auxv , dim);				//x = x + alpha * p;

		rnew = inprodvr(r,r,dim);				//rsnew = r' * r;

		printf("\r Iter = %d ; |r_CG| = %.0e ",(*iCG),sqrt(rnew));
		fflush(stdout);

		if( (sqrt(rnew) < tol) || (sqrt(rnew) > 1e2) || *iCG > ITERATION_TOL)					//if sqrt(rsnew) < tol
			break;								//break

		cprodvr(auxv , rnew/rold , p ,dim);
		sumvr(p , r , auxv , dim);				//p = r + (rsnew / rsold) * p;

		rold = rnew;							//rsold = rsnew;
	}

	free(r);
	free(p);
	free(auxv);

	*r_CG = rnew;

	return((double)(clock()-dtime)/CLOCKS_PER_SEC);
}

void rgetOrthogonalSpace(double * c, double * b, LatticeSU2 * lattice, double cg_tol, int * iCG, double * r_CG){
	double * a = malloc(sizeof(double)*colorV);

	//getting the null space out of M
	rfpapply(a , lattice , b , colorV);
	//solve CG
	copyvr(c,a,colorV);
	//double pmin[4] = {0e0,0e0,0e0,1e0/N};
	//rsetplanewave(c, pmin ,0);
	//setonevc(c,colorV);
	rfpcgmr(c, lattice , a, colorV, cg_tol, iCG, r_CG);

	free(a);
}

int rsmallest_eigen_cg(LatticeSU2 * lattice, double * lambda1, double * eigen_out
	 	, double eigen_tol , double * vguess, double * r_cg, int * iCG, double cg_tol){
	int i;
	clock_t dtime = clock();
	double norm;
	double eigen_v = 0e0;
	double eigen_vp = 0e0;
	double * X = malloc(sizeof(double)*colorV);
	double * Y = malloc(sizeof(double)*colorV);
	double * Xp = malloc(sizeof(double)*colorV);
	*r_cg = 0e0;
	//double pmin[4] = {0e0,0e0,0e0,1e0/N};
	//rsetplanewave(X, pmin ,2);
	//rsetplanewaveallcolors(X, pmin);
	//setonevr(X,colorV);
	//setrandomvr(X, colorV, _seed);

	copyvr(X,vguess,colorV);											//vguess is some guess to the eigenvector with the smallest eigenvalue

	//Power method
	i=0;

	//rgetOrthogonalSpace(X,X,lattice);
	printf("\nPM: Iter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/eigen_vp,norm);
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
		copyvr(Xp,X,colorV);

		rfpcgmr(Xp , lattice , X , colorV, cg_tol, iCG, r_cg);
		if(*r_cg > 1e2 || *iCG>ITERATION_TOL){
			printf("\n-----> r_CG DIVERGED in PM\n");
			fflush(stdout);
			//getchar();
			break;
		}

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodvr(X,Xp,colorV)/inprodvr(X,X,colorV);
		copyvr(X,Xp,colorV);
		reunitvr(X,colorV);

		i++;
		eigen_vp = (eigen_vp);
		rfpapply(Xp,lattice,X,colorV);
		cprodvr(Y,-1e0/eigen_vp,X,colorV);
		sumvr(Xp,Xp,Y,colorV);
		norm = normvr(Xp,colorV);
		printf("\nPM: Iter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/eigen_vp,norm);
		fflush(stdout);

		cprodvr(X, getSignal(min(X,colorV)), X, colorV);		//Attilio's trick

	}while( norm > eigen_tol);	//this way the eigenvalue equation is directly verified, but is slower

	//	printf("\n ev = %lf | evp = %lf | abs() = %lf | tol = %lf \n",1.0/eigen_v,1.0/eigen_vp,fabs(1.0/eigen_v - 1.0/eigen_vp),eigen_tol);
	//}while( fabs(1.0/eigen_v - 1.0/eigen_vp) > 1e-7 ); 	//this way is faster but the tolerance must be smaller

	copyvr(eigen_out, X, colorV);

	free(X);
	free(Xp);
	free(Y);
	*lambda1 = 1e0/eigen_vp;
	return( i );
}

//CG and PM with complex vector
double fpapply(double complex * Mvout , LatticeSU2 * lattice , double complex * vin , int dim){
	//dim should be = N_colors * N^4
	//the color index = a is fixed
	//ATTENTION: a is 1,2,3 or 0,1,2?
	int a,b,c,t,x,y,z,mi;
	int ami[4];
	double complex * Mvouttemp;
	clock_t dtime = clock();
	Mvouttemp = malloc(sizeof(double complex)*dim);

	setzerovc(Mvouttemp,dim);

	for(a=0;a<3;a++){
		for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
		for(y=0;y<N;y++){
		for(z=0;z<N;z++){
			for(mi=0;mi<4;mi++){

				setzeroqvector(ami);
				ami[mi]=1;

				*getelvc(Mvouttemp,a,t,x,y,z) +=

				*getelv(getU(lattice,t,x,y,z,mi),0)
				*(
				*getelvc(vin,a,t,x,y,z)
				-*getelvc(vin,a,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
				)

				+

				*getelv(getU(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi),0)
				*(
				*getelvc(vin,a,t,x,y,z)
				-*getelvc(vin,a,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
				);

				//sum in b,c
				for(b=0;b<3;b++){
					for(c=0;c<3;c++){
						if((a!=b) && (b!=c)){
							*getelvc(Mvouttemp,a,t,x,y,z) +=
							aseps(a,b,c)*
							(

							(*getelv(getU(lattice,t,x,y,z,mi),b+1))*
							(
								*getelvc(vin,c,getStepT(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]))
								-
								0//*getelvc(vin,c,t,x,y,z)
							)

							-

							(*getelv(getU(lattice,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]),mi),b+1))*
							(
								*getelvc(vin,c,getStepT(t,-ami[0]),getStep(x,-ami[1]),getStep(y,-ami[2]),getStep(z,-ami[3]))
								-
								0//*getelvc(vin,c,t,x,y,z)
							)
							);
						}
					}
				}
			}
		}}}}
	}

	copyvc(Mvout,Mvouttemp,dim);
	free(Mvouttemp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

double fpcgmr(double complex * x , LatticeSU2 * lattice , double complex * b,  int dim){
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a hermitian (dim x dim) matrix (real and symmetric in this case)
	//A is the FP matrix defined from the lattice links applied to the source b
	//in the usual context dim = N_colors * N^4 ==> the vector will have a super-index (a,x)

	double tol = 1E-10;
	int i,icont = 0;
	double complex alpha, rnew,rold;;
	double complex * auxv = malloc(sizeof(double complex)*dim);
	double complex * r = malloc(sizeof(double complex)*dim);
	double complex * p = malloc(sizeof(double complex)*dim);
	clock_t dtime = clock();

	//chute inicial passado em x
	//definições iniciais:

	fpapply(auxv , lattice , x  ,dim);		//mvprodr(auxv , A , x , dim);

	//fpapply(auxv , lattice , auxv  ,dim);		//mvprodr(auxv , A , x , dim);
	//fpapply(b , lattice , b  ,dim);		//mvprodr(auxv , A , x , dim);

	cprodvc(auxv , -1e0 , auxv , dim);
	sumvc(r , b , auxv , dim);				//r <- b - A.x

	copyvc(p , r,dim);						//p <- r

	rold = inprodvc(r,r,dim);

	printf("\n|r_CG| = %.0e ",sqrt(rnew));
	fflush(stdout);

	while(1){														//for i = 1:length(b)
		icont++;

		fpapply(auxv,lattice,p,dim);			//Ap = A * p;

		alpha = rold/inprodvc(p,auxv,dim);		//alpha = rsold / (p' * Ap);

		cprodvc(auxv , -alpha , auxv, dim);
		sumvc(r , r , auxv , dim);				//r = r - alpha * Ap;

		cprodvc(auxv , alpha , p ,dim);
		sumvc(x , x , auxv , dim);				//x = x + alpha * p;

		rnew = inprodvc(r,r,dim);				//rsnew = r' * r;

		printf("\r|r_CG| = %.0e ",sqrt(rnew));
		fflush(stdout);

		if(sqrt(rnew) < tol)					//if sqrt(rsnew) < tol
			break;								//break

		cprodvc(auxv , rnew/rold , p ,dim);
		sumvc(p , r , auxv , dim);				//p = r + (rsnew / rsold) * p;

		rold = rnew;							//rsold = rsnew;
	}

	free(r);
	free(p);
	free(auxv);

	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

void getOrthogonalSpace(double complex * c, double complex * b, LatticeSU2 * lattice){
	double complex * a = malloc(sizeof(double complex)*colorV);
	//getting the null space out of M
	fpapply(a , lattice , b , colorV);
	//solve CG
	copyvc(c,a,colorV);
	//setonevc(c,colorV);
	fpcgmr(c, lattice , a, colorV);

	free(a);
}

double complex smallest_eigen_cg(LatticeSU2 * lattice, double * lambda1, double complex * eigen_out){
	int i;
	double norm;
	double eigen_tol = 1e-5;
	clock_t dtime = clock();
	double complex eigen_v = 0e0 + I*0e0;
	double complex eigen_vp = 0e0 + I*0e0;
	double complex * X = malloc(sizeof(double complex)*colorV);
	double complex * Y = malloc(sizeof(double complex)*colorV);
	double complex * Xp = malloc(sizeof(double complex)*colorV);
	double pmin[4] = {0e0,0e0,0e0,1e0/N};
	setplanewaveallcolors(X, pmin);

	//Power method
	i=0;
	printf("\nIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);
	do{
		getOrthogonalSpace(X,X,lattice);
		copyvc(Xp,X,colorV);
		fpcgmr(Xp , lattice , X , colorV);

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodvc(X,Xp,colorV)/inprodvc(X,X,colorV);
		copyvc(X,Xp,colorV);
		reunitvc(X,colorV);

		i++;

		fpapply(Xp,lattice,X,colorV);
		cprodvc(Y,-1e0/eigen_vp,X,colorV);
		sumvc(Xp, Xp,Y,colorV);

		norm = normvc(Xp,colorV);

		printf("\rIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);

	}while( norm > eigen_tol );
	copyvc(eigen_out, X, colorV);

	free(X);
	free(Xp);
	free(Y);
	*lambda1 = 1e0/creal(eigen_vp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

//biCGstab with complex vectors
void fpbicgstabmr(double complex * x , LatticeSU2 * lattice , double complex * b,  int dim, double cg_tol, double * r_cg, int * icg){
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
	double complex * aux1 = malloc(sizeof(double complex)*dim);
	double complex * aux2 = malloc(sizeof(double complex)*dim);
	double complex * rtilde = malloc(sizeof(double complex)*dim);
	double complex * r = malloc(sizeof(double complex)*dim);
	double complex * v = malloc(sizeof(double complex)*dim);
	double complex * s = malloc(sizeof(double complex)*dim);
	double complex * t = malloc(sizeof(double complex)*dim);
	double complex * p = malloc(sizeof(double complex)*dim);
	clock_t dtime = clock();

	//chute inicial passado em x
	//definições iniciais:

	fpapply(aux1 , lattice , x  ,dim);		//aux1 <- M.x

	cprodvc(aux1 , -1e0 , aux1 , dim);
	sumvc(r , b , aux1 , dim);				//r <- b - A.x

	copyvc(rtilde, r,dim);						//rtilde <- r

	printf("\n iter=%d |s_biCGstab| = %e", *icg, *r_cg);
	fflush(stdout);
	while(1){//for(i=0;i<dim;i++){				//for i = 1:length(b)
		(*icg)++;

		rho1 = inprodvc(rtilde, r , dim);

		if( cabs(rho1) == 0){
			printf("\nrho=0 : BiCGStab failed.");
			break;
		}

		if(*icg==1){
			copyvc(p,r,dim);
		}
		else{
			beta = alpha*rho1/(omega*rho2);

			cprodvc(aux1,-omega,v,dim);
			sumvc(aux1,aux1,p,dim);
			cprodvc(aux1,beta,aux1,dim);
			sumvc(p,r,aux1,dim);			//p <- r + beta.(p-omega.v)
		}

		fpapply(v,lattice,p,dim);			//v <- A.p

		alpha = rho1/inprodvc(rtilde,v,dim);	//alpha <- rho1/(rtilde,v)

		cprodvc(aux1,-alpha,v,dim);
		sumvc(s,r,aux1,dim);

		*r_cg = normvc(s,dim);

		printf("\r iter=%d |s_biCGstab| = %e", *icg, *r_cg);
		fflush(stdout);

		if( *r_cg < tol || *r_cg > 1e2 || *icg>ITERATION_TOL){
			cprodvc(aux1,alpha,p,dim);
			sumvc(x,x,aux1,dim);			//x <- x + alpha.p
			break;
		}

		fpapply(t,lattice,s,dim);			//t <- A.s

		omega = inprodvc(t,s,dim)/inprodvc(t,t,dim);	//omega <- (t,s)/(t,t)

		cprodvc(aux1,-omega,t,dim);
		sumvc(r,s,aux1,dim);				//r <- s - omega.t

		cprodvc(aux1,alpha,p,dim);
		cprodvc(aux2,omega,s,dim);
		sumvc(aux1,aux1,aux2,dim);
		sumvc(x,x,aux1,dim);				//x <- x + alpha.p + oemga.s

		rho2 = rho1;
	}

	free(r);
	free(rtilde);
	free(v);
	free(s);
	free(t);
	free(p);
	free(aux1);
	free(aux2);

	//return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

void bigetOrthogonalSpace(double complex * c, double complex * b, LatticeSU2 * lattice, double cg_tol, double * r_cg, int * icg){
	double complex * a = malloc(sizeof(double complex)*colorV);
	//getting the null space out of M
	fpapply(a , lattice , b , colorV);
	//solve CG
	copyvc(c,a,colorV);
	//setonevc(c,colorV);
	fpbicgstabmr(c, lattice , a, colorV,cg_tol,r_cg,icg);

	free(a);
}

double smallest_eigen_bicgstab(LatticeSU2 * lattice, double * lambda1, double * eigen_out, double eigen_tol, double * r_cg, int * icg, double cg_tol){
	int i,t,x,y,z,a;
	clock_t dtime = clock();
	double norm;
	double complex eigen_v = 0e0 + I*0e0;
	double complex eigen_vp = 0e0 + I*0e0;
	double complex * X = malloc(sizeof(double complex)*colorV);
	double complex * Y = malloc(sizeof(double complex)*colorV);
	double complex * Xp = malloc(sizeof(double complex)*colorV);
	double pmin[4] = {0e0,0e0,0e0,1e0/N};

	//converting to complex vector_guess, i.e. to X
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					for(a=0;a<3;a++){
						*getelvc(X,a,t,x,y,z) = *getelvr(eigen_out,a,t,x,y,z);
					}
				}
			}
		}
	}


	//Power method
	i=0;
	printf("\nIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);
	fflush(stdout);
	do{
		bigetOrthogonalSpace(X,X,lattice,cg_tol,r_cg,icg);
		if(*r_cg > 1e2 || *icg>ITERATION_TOL){
			printf("\n-----> biCG DIVERGED in PM (orthogonalization)! ABORTING CG AND TRYING NEW GUESS\n");
			fflush(stdout);
			//getchar();
			break;
		}
		copyvc(Xp,X,colorV);
		fpbicgstabmr(Xp , lattice , X , colorV,cg_tol,r_cg,icg);
		if(*r_cg > 1e2 || *icg>ITERATION_TOL){
			printf("\n-----> biCG DIVERGED in PM! ABORTING CG AND TRYING NEW GUESS\n");
			fflush(stdout);
			//getchar();
			break;
		}

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodvc(X,Xp,colorV)/inprodvc(X,X,colorV);
		copyvc(X,Xp,colorV);
		reunitvc(X,colorV);

		i++;

		fpapply(Xp,lattice,X,colorV);
		cprodvc(Y,-1e0/eigen_vp,X,colorV);
		sumvc(Xp, Xp,Y,colorV);

		norm = normvc(Xp,colorV);

		printf("\rIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);
		fflush(stdout);

	}while( (norm > eigen_tol)  );
	//copyvc(eigen_out, X, colorV);

	//converting back to real output vector
	for(t=0;t<Nt;t++){
		for(x=0;x<N;x++){
			for(y=0;y<N;y++){
				for(z=0;z<N;z++){
					for(a=0;a<3;a++){
						*getelvr(eigen_out,a,t,x,y,z) = creal(*getelvc(X,a,t,x,y,z));
					}
				}
			}
		}
	}
	reunitvr(eigen_out,colorV);

	free(X);
	free(Xp);
	free(Y);
	*lambda1 = 1e0/creal(eigen_vp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

//R functional computation
double  FsecondDerivative(LatticeSU2 * lattice, double complex * eigen_vector_out){
	double lambda;
	smallest_eigen_cg(lattice,&lambda,eigen_vector_out);
	return lambda/(4e0*totalV);
}

double rFThirdFourthDerivative(LatticeSU2 * lattice, double * ev_in, double * third, double * fourth){
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
						setquadv(emi,mi);
						for(a=0;a<3;a++){
							for(b=0;b<3;b++){

								*third -=
								0.75*( pow(*getelvr(ev_in,a,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3])),2)
									    - pow(*getelvr(ev_in,a,t,x,y,z),2) )
									*( *getelvr(ev_in,b,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3]))
									    + *getelvr(ev_in,b,t,x,y,z) )
									*( getU(lattice,t,x,y,z,mi)[b+1] );

								*third -=
								 (pow(*getelvr(ev_in,a,t,x,y,z),2) )
									*( *getelvr(ev_in,b,t,x,y,z) )
									*( getU(lattice,t,x,y,z,mi)[b+1]
										- getU(lattice,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]),mi)[b+1] );


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
								*getU(lattice,t,x,y,z,mi)[0];

								*fourth -= pow(*getelvr(ev_in,a,t,x,y,z),2)
											*( *getelvr(ev_in,b,t,x,y,z) )*
									(
										getU(lattice,t,x,y,z,mi)[0]*(
											*getelvr(ev_in,b,t,x,y,z) -  *getelvr(ev_in,b,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3]))
										)

										+

										getU(lattice,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]),mi)[0]*(
											*getelvr(ev_in,b,t,x,y,z) - *getelvr(ev_in,b,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]))
										)
									);

								for(d=0;d<3;d++){
									for(e=0;e<3;e++){
										if((b!=d) && (d!=e) && (b!=e)){
											*fourth -= pow(*getelvr(ev_in,a,t,x,y,z),2)*( *getelvr(ev_in,b,t,x,y,z) )*
												(
													aseps(b,d,e)*(
														getU(lattice,t,x,y,z,mi)[d+1]*( *getelvr(ev_in,e,getStepT(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3])) )
														-
														getU(lattice,getStepT(t,-emi[0]),getStep(x,-emi[1]),getStep(y,-emi[2]),getStep(z,-emi[3]),mi)[d+1]*(
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
double complex largest_eigen_cg(LatticeSU2 * lattice, double * lambda1, double complex * eigen_out){
	int i;
	double eigen_tol = 1e-8;
	double complex eigen_v = 0e0 + I*0e0;
	double complex eigen_vp = 0e0 + I*0e0;
	double complex * X = malloc(sizeof(double complex)*colorV);
	double complex * Y = malloc(sizeof(double complex)*colorV);
	double complex * Xp = malloc(sizeof(double complex)*colorV);
	double pmin[4] = {0e0,0e0,0e0,1e0/N};
	setplanewaveallcolors(X, pmin);
	//setonevc(X,colorV);
	clock_t dtime = clock();

	double norm;

	//Power method
	i=0;
	//printf("\nIter=%d | lambda_1=%e |  lambda1_norm=%e",i,1e0/creal(eigen_vp),norm);
	do{

		//CG
		//setplanewaveallcolors(Xp, pmin);
		//setonevc(Xp,colorV);
		//setzerovc(Xp,colorV);
		fpapply(Xp , lattice , X , colorV);

		//PM
		eigen_v = eigen_vp;
		eigen_vp = inprodvc(X,Xp,colorV)/inprodvc(X,X,colorV);
		copyvc(X,Xp,colorV);
		reunitvc(X,colorV);

		i++;

		fpapply(Xp,lattice,X,colorV);
		cprodvc(Y,-eigen_vp,X,colorV);
		sumvc(Xp, Xp,Y,colorV);

		norm = normvc(Xp,colorV);

		//printf("\rIter=%d | lambda_1=%e |  lambda1_norm=%e",i,creal(eigen_vp),norm);

	}while( norm > eigen_tol );
	copyvc(eigen_out, X, colorV);

	free(X);
	free(Xp);
	free(Y);
	*lambda1 = 1e0/creal(eigen_vp);
	return( (double)(clock()-dtime)/CLOCKS_PER_SEC );
}

void divergentConfigProcedure(LatticeSU2 * escaled_lattice, double new_delta_scale, double * tau, int * divergenConfigN, double r_CG, int iCG
																			,double * lambda1, double * ev_guess, double eigen_tol, double cg_tol){
	//in case of patological configuration
	double pmin[4] = {0,0,0,1.0/N};
	(*divergenConfigN)++;

	while(r_CG > 1e2 || iCG > ITERATION_TOL){
		if(r_CG > 1e2 || iCG > ITERATION_TOL){
			//if none of the above worked, try scaling back and using a step midways and repeat all
			*tau /= new_delta_scale;
			reescaleGaugeField(escaled_lattice, escaled_lattice, 1.0/new_delta_scale);		//cancelling the last scaling
			new_delta_scale = (1.0+new_delta_scale)/2.0;
			*tau *= new_delta_scale;
			reescaleGaugeField(escaled_lattice, escaled_lattice, new_delta_scale);		//and rescaling midways
			printf("\nTRYING MIDWAY STEP=%e",new_delta_scale);
		}
		if(r_CG > 1e2 || iCG > ITERATION_TOL){
			//setting a random start and trying again
			printf("\nTRYING CG WITH RANDOM GUESS");
			setrandomvr(ev_guess,colorV);
			rsmallest_eigen_cg(escaled_lattice, lambda1, ev_guess, eigen_tol, ev_guess , &r_CG, &iCG, cg_tol);
		}
		if(r_CG > 1e2 || iCG > ITERATION_TOL){
			//if the last try didnt work, try using bicg (with plane wave first)
			printf("\nTRYING biCG WITH PW GUESS");
			rsetplanewaveallcolors(ev_guess, pmin);
			smallest_eigen_bicgstab(escaled_lattice, lambda1, ev_guess, eigen_tol, &r_CG, &iCG, cg_tol);
		}
		if(r_CG > 1e2 || iCG > ITERATION_TOL){
			//if the last try didnt work, try using bicg (with random vector second)
			printf("\nTRYING biCG WITH RANDOM GUESS");
			setrandomvr(ev_guess,colorV);
			smallest_eigen_bicgstab(escaled_lattice, lambda1, ev_guess, eigen_tol, &r_CG, &iCG, cg_tol);
		}
		//save_lattice( target_lattice , "./bin/divlattice.bin");
		fflush(stdout);
	}
}

double horizonWalk(LatticeSU2 * target_lattice, double * nHorizonOut, double * rhoHorizonOut, double * ev_guess
		, double * lambda1Out, double * rAfterOut, double * rBeforeOut, double * thirdOut, double * thirdAbsOut, double * fourthOut
		, double * rOut, double * pwProjOut, double eigen_tol, int * divergentConfig, double cg_tol){

	//CG INVERSION AND DATA COLLECTION OF THE FIRST STEP
	double delta_scale=1e0, tau=1e0;
	double r_CG = 0;
	int iCG=0;

	rsmallest_eigen_cg(target_lattice, lambda1Out, ev_guess, eigen_tol, ev_guess, &r_CG, &iCG, cg_tol);
	if(r_CG > 1e2 || iCG > ITERATION_TOL){
		printf("\n--> DIVERGED ON START! r_CG=%e | iCG=%d | LAST: lambda1=%e",r_CG,iCG,*lambda1Out);
		divergentConfigProcedure(target_lattice, delta_scale, &tau, divergentConfig, r_CG, iCG, lambda1Out, ev_guess, eigen_tol, cg_tol);
	}
	rFThirdFourthDerivative(target_lattice, ev_guess, thirdOut, fourthOut);
	*thirdAbsOut = fabs(*thirdOut);
	*rOut = 4*totalV*(*thirdOut)*(*thirdOut)/((*fourthOut)*(*lambda1Out));
	*pwProjOut = compareVectorToPW(ev_guess , 1e0/N);
	printf("\n Initial: %e %e %e %e\n",*lambda1Out,*thirdOut,*fourthOut,*rOut);	fflush(stdout);

	printf("Initial data save. Starting walk...\n");

	// double * escaled_lattice = malloc(sizeof(double)*dimLattice);
	LatticeSU2 * escaled_lattice = newLatticeSU2(target_lattice->N , target_lattice->Nt);

	double lambda1, third, fourth, r;
	*nHorizonOut = 0;
	do{
		//gauge field reescaling selection
		if(lambda1 > 5e-3)
			delta_scale = 1.001;
		else if(5e-3 >= lambda1 && lambda1 > 5e-4)
			delta_scale = 1.0005;
		else
			delta_scale = 1.0001;

		tau *= delta_scale;

		printf("\n nHorizon=%d ; Reescaling Gauge field by tau=%lf\n",(int)*nHorizonOut,tau);
		reescaleGaugeField(escaled_lattice, target_lattice, tau);
		rsmallest_eigen_cg(escaled_lattice, &lambda1, ev_guess, eigen_tol, ev_guess , &r_CG, &iCG, cg_tol);
		if(r_CG > 1e2 || iCG > ITERATION_TOL){
			printf("\n--> DIVERGED! r_CG=%e | iCG=%d | LAST r=%e , E'''=%e",r_CG,iCG,r,third);
			divergentConfigProcedure(escaled_lattice, delta_scale, &tau, divergentConfig, r_CG, iCG, &lambda1, ev_guess, eigen_tol, cg_tol);
		}
		rFThirdFourthDerivative(escaled_lattice, ev_guess, &third, &fourth);
		*rBeforeOut=r;
		r = 4.0*totalV*third*third/(fourth*lambda1);
		(*nHorizonOut)++;

		printf("\nN=%d | Horizon step=%d %e %e %e %e %e %e\n",N,(int)*nHorizonOut,lambda1,third,fourth,r,delta_scale,tau);
		fflush(stdout);
	}while(lambda1 > 0e0);
	*rAfterOut = r;
	*rhoHorizonOut = 1e0/((1e0/delta_scale + 1e0)*tau/2e0);
	printf(" -- rho = %lf",*rhoHorizonOut);

	free(escaled_lattice->U);
	return(r_CG);
}

double verifyOrthogonalization(LatticeSU2 * lattice, double * vin){
	double * aux 		= malloc(sizeof(double)*N*N*N*N*3);
	rfpapply(aux, lattice, vin, colorV);

	double * constant_vector 		= malloc(sizeof(double)*N*N*N*N*3);
	setonevr(constant_vector,colorV);

	double result = inprodvr(constant_vector,aux,colorV);

	free(constant_vector);
	free(aux);

	return(result);
}

