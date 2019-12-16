#include "algebra.h"
#include "global.h"


int signal(double number){
	//returns +1 if element is (zero) positive and -1 else
	if(number < 0)	return -1;
	return 1;
}

////////////////////////////// SU2 MATRIX/VECTOR ALGEBRA

double tracev( double * pvector ){
	int i;
	return 2*(*getelv(pvector, 0));
}

void sumv( double * pout , double * B , double * C ){
	double  A[4];
	int i,j;

	for(i=0 ; i<4 ; i++){
		A[i] = *getelv(B,i) + *getelv(C,i);
	}

	for(i=0 ; i<4 ; i++){
			*getelv(pout,i) = A[i];
	}
}

void mmprodv( double * pout , double  * B, double  * C ){
	double A[4];
	int i;

	A[0] =   (*(getelv(B,0)) )*( *(getelv(C,0)))-
			 (*(getelv(B,3)) )*( *(getelv(C,3)))-
			 (*(getelv(B,2)) )*( *(getelv(C,2)))-
			 (*(getelv(B,1)) )*( *(getelv(C,1)));

	A[1] =  -(*(getelv(C,3)) )*( *(getelv(B,2)))+
			 (*(getelv(C,0)) )*( *(getelv(B,1)))+
			 (*(getelv(C,1)) )*( *(getelv(B,0)))+
			 (*(getelv(C,2)) )*( *(getelv(B,3)));

	A[2] =   (*(getelv(C,2)) )*( *(getelv(B,0)))-
			 (*(getelv(C,1)) )*( *(getelv(B,3)))+
			 (*(getelv(C,0)) )*( *(getelv(B,2)))+
			 (*(getelv(C,3)) )*( *(getelv(B,1)));

	A[3] =   (*(getelv(C,3)) )*( *(getelv(B,0)))+
			 (*(getelv(C,0)) )*( *(getelv(B,3)))+
			 (*(getelv(C,1)) )*( *(getelv(B,2)))-
			 (*(getelv(C,2)) )*( *(getelv(B,1)));

	for(i=0 ; i<4 ; i++){
			*(getelv(pout,i)) = A[i];
	}
}

//need revising
void cmprodv( double * pout , double k, double * C ){
	//double A[4];
	//double a = creal(k);
	//double b = cimag(k);
	int i;

	//A[0] = a*(*getelv(C,0)) - b*(*getelv(C,3));
	//A[1] = a*(*getelv(C,1)) - b*(*getelv(C,2));
	//A[2] = a*(*getelv(C,2)) + b*(*getelv(C,1));
	//A[3] = a*(*getelv(C,3)) + b*(*getelv(C,0));

	for(i=0 ; i<4 ; i++){
		*(getelv(pout,i)) = k*(*(getelv(C,i)));
	}
}
//

void hermcv(double * Adagger , double * A){
	int i;
	*getelv(Adagger,0) = *getelv(A,0);
	for(i=1 ; i<4 ; i++){
		*getelv(Adagger,i) = - *getelv(A,i) ;
	}
}

double detv(double * A){
	return (*getelv(A,0)*(*getelv(A,0))
		  + *getelv(A,1)*(*getelv(A,1))
		  + *getelv(A,2)*(*getelv(A,2))
		  + *getelv(A,3)*(*getelv(A,3)));
}

void setidv(double * vector){
	int i;
	*getelv(vector , 0) = 1.0;
	for (i=1 ; i<4 ; i++){
		*getelv(vector , i) = 0.0;
	}
	//printVector(vector);
}

void matrixexpv(double * vout, double * vin, int order){
	int i,j;
	double fat= 1E0;
	double * outtemp = malloc(sizeof(double)*4);
	double *aux = malloc(sizeof(double)*4);
	double *A = malloc(sizeof(double)*4);
	setidv(outtemp);
	setidv(A);

	for(i=1;i <= order;i++){
		mmprodv(A,vin,A);
		cmprodv(aux,1/fat,A);


		//printf("\n%d , %lf\n",i,1/fat);
		//printm(aux);
		//printf("\n----\n");
		sumv(outtemp,outtemp,aux);
		//printm(outtemp);
		fat *= i+1E0;
	}
	copyv(vout,outtemp);
	free(outtemp);
	free(aux);
	free(A);
}

void setzerov(double * vector){
	int i;
	for (i=0 ; i<4 ; i++){
		*getelv(vector , i) = 0.0;
	}
}

void su2sroot(double * out, double * in){
	int i;
	double * temp = malloc(sizeof(double)*4);

	reunitv(in);

	temp[0] = cos(acos(in[0])/2e0);
	for(i=1;i<4;i++){
		temp[i] = in[i]/(2e0*temp[0]);
	}
	copyv(out,temp);

	free(temp);
}

////////////////////////////// GENERAL REAL MATRIX ALGEBRA

double tracer( double * pmatrix , int dim){
	int i;
	double trace = 0E0;
	for(i=0;i<dim;i++){
		trace += *getelmr(pmatrix, i , i , dim);
	}
	return (trace);
}

void setidmr(double * M , int dim){
	int i,j;

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmr(M,i,j,dim) = 1E0*(i==j) ;
		}
	}
}

void setzeromr(double * M , int dim){
	int i,j;

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmr(M,i,j,dim) = 0E0 ;
		}
	}
}

void sumr( double * pout , double * B , double * C , int dim){
	int i,j;

	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			*getelmr(pout,i,j,dim) = *getelmr(B,i,j,dim) + *getelmr(C,i,j,dim);
		}
	}
}

void mmprodr( double * mout , double * B , double * C , int dim){
	double * Atemp = malloc(sizeof(double)*dim*dim);
	int i,j,k;

	setzeromr(Atemp,dim);
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			for(k=0 ; k<dim ; k++){
				*getelmr(Atemp,i,j,dim) +=  *(getelmr(B,i,k,dim)) * (*getelmr(C,k,j,dim));
			}
		}
	}

	copymr(mout,Atemp,dim);

	free(Atemp);
}

void transpmr(double * Atransposed , double * A , int dim){
	//computes the transpose of the real matrix A
	int i,j;
	double * Atemp = malloc(sizeof(double)*dim*dim);

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmr(Atemp,i,j,dim) = *getelmr(A,j,i,dim);
		}
	}

	copymr(Atransposed , Atemp , dim);
	free(Atemp);
}

void setzerovr(double * vr , int dim){
	//dim must be the dimension of the vector vc of size (double complex)*dim
	int i;
	for(i=0;i<dim;i++){
		vr[i] = 0E0;
	}
}

double inprodvr(double * a , double * b , int dim){
	//compute the usual inner product defined in a vetorial space of real variables
	//returns (a,b) = a . b
	int i;
	double sum = 0E0;
	for(i=0;i<dim;i++){
		sum += a[i] * b[i];
	}
	return(sum);
}

double normvr(double * a , int dim){
	//returns norm defined by inprov : |a| = sqrt(a,a)
	return( sqrt(inprodvr(a,a,dim)) );
}

void cprodvr(double * vout , double k , double * v , int dim){
	int i;

	for(i=0;i<dim;i++){
		vout[i] = k*v[i];
	}
}

void sumvr(double * vout , double * v1 , double * v2 , int dim){
	int i;

	for(i=0;i<dim;i++){
		vout[i] = v1[i] + v2[i];
	}
}

void reunitvr(double * vector, int dim){
	int i;
	double norm;
	norm = normvr(vector, dim);
	for(i=0;i<dim;i++){
		vector[i] /= norm;
	}
}

void mvrprod(double * vout , double * B , double * vcolumn , int dim){
	//computes vout  <- B . Vcolumn
	//B is a real (dim x dim) matrix
	//vcolumn is a real (dim x 1) column vector
	//vout is a real (dim x 1) column vector
	int i,j;
	double * vtemp = malloc(sizeof(double)*dim);

	setzerovr(vtemp,dim);
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			vtemp[i] +=  *getelmr(B,i,j,dim)*vcolumn[j];
		}
	}

	copyvr(vout,vtemp,dim);
	free(vtemp);
}

////////////////////////////// GENERAL COMPLEX MATRIX ALGEBRA

void setzerovc(double complex * vc , int dim){
	//dim must be the dimension of the vector vc of size (double complex)*dim
	int i;
	for(i=0;i<dim;i++){
		vc[i] = 0E0;
	}
}

double complex tracec( double complex * pmatrix , int dim){
	int i;
	double complex trace = 0E0;
	for(i=0;i<dim;i++){
		trace += *getelmc(pmatrix, i , i , dim);
	}
	return (trace);
}

void setidmc(double complex * M , int dim){
	int i,j;

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmc(M,i,j,dim) = 1E0*(i==j) + 0E0*I ;
		}
	}
}

void setzeromc(double complex * M , int dim){
	int i,j;

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmc(M,i,j,dim) = 0E0 + 0E0*I ;
		}
	}
}

void sumc( double complex * mout , double complex * B , double complex * C , int dim){
	int i,j;

	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			*getelmc(mout,i,j,dim) = *getelmc(B,i,j,dim) + *getelmc(C,i,j,dim);
		}
	}
}

void mmprodc( double complex * Mout , double complex * B , double complex * C , int dim){
	int i,j,k;
	double complex * Mtemp = malloc(sizeof(double complex)*dim*dim);

	setzeromc(Mtemp,dim);
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			for(k=0 ; k<dim ; k++){
				*getelmc(Mtemp,i,j,dim) +=  *(getelmc(B,i,k,dim)) * (*getelmc(C,k,j,dim));
			}
		}
	}
	copymc(Mout,Mtemp,dim);

	free(Mtemp);
}

void cmprodc(double complex * mout , double complex k , double complex * C , int dim){
	int i,j;

	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			*getelmc(mout,i,j,dim) = k*( *getelmc(C,i,j,dim));
		}
	}
}

void mvprodr(double complex * vout , double * B , double complex * vcolumn , int dim){
	//computes vout  <- B . Vcolumn
	//B is a real (dim x dim) matrix
	//vcolumn is a complex (dim x 1) column vector
	//vout is a complex (dim x 1) column vector
	int i,j;
	double complex * vtemp = malloc(sizeof(double complex)*dim);

	setzerovc(vtemp,dim);
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			vtemp[i] +=  *getelmr(B,i,j,dim)*vcolumn[j];
		}
	}

	copyvc(vout,vtemp,dim);
	free(vtemp);
}

void vmprodr(double complex * vout , double complex* vline , double * B , int dim){
	//computes vout  <- vline . B
	//B is a real (dim x dim) matrix and
	//vline is a complex (1 x dim) line vector
	//vout is a complex (1 x dim) line vector
	int i,j;
	double complex * vtemp = malloc(sizeof(double complex)*dim);

	setzerovc(vtemp,dim);
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			vtemp[i] +=  *getelmr(B,i,j,dim)*vline[j];
		}
	}

	copyvc(vout,vtemp,dim);
	free(vtemp);
}

void sumvc(double complex * vout , double complex * v1 , double complex * v2 , int dim){
	int i;

	for(i=0;i<dim;i++){
		vout[i] = v1[i] + v2[i];
	}
}

void cprodvc(double complex * vout , double complex k , double complex * v , int dim){
	int i;

	for(i=0;i<dim;i++){
		vout[i] = k*v[i];
	}
}

void hermcvc(double complex * Vdagger , double complex * V , int dim){
	//computes the hermitian conjugate of the complex vector V
	int i;
	double complex * Vtemp = malloc(sizeof(double complex)*dim);

	for(i=0;i<dim;i++){
		Vtemp[i] = conj( V[i] );
	}

	copyvc(Vdagger , Vtemp , dim);
	free(Vtemp);
}

void hermcmc(double complex * Mdagger , double complex * M , int dim){
	//computes the hermitian conjugate of the complex matrix M
	int i,j;
	double complex * Mtemp = malloc(sizeof(double complex)*dim*dim);

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmc(Mtemp,i,j,dim) = conj( *getelmc(M,j,i,dim) );
		}
	}

	copymc(Mdagger , Mtemp , dim);
	free(Mtemp);
}

double complex inprodvc(double complex * a , double complex * b , int dim){
	//compute the usual inner product defined in a vetorial space of complex variables
	//returns (a,b) = adagger . b
	int i;
	double complex sum = 0E0 + 0E0*I;
	double complex * adagger = malloc(sizeof(double complex)*dim);
	hermcvc(adagger , a , dim);
	//printf("\n\n adagger:");
	//printvc(b,dim);
	for(i=0;i<dim;i++){
		sum += adagger[i] * b[i];
	}
	//printf("\n\n inner product:");
	//printc(sum);
	free(adagger);
	return(sum);
}

double normvc(double complex * a , int dim){
	//returns norm defined by inprovc : |a| = sqrt(a,a)
	return( creal(sqrt(inprodvc(a,a,dim))) );
}

void reunitvc(double complex * vector, int dim){
	int i;
	double norm;
	norm = normvc(vector, dim);
	for(i=0;i<dim;i++){
		vector[i] /= norm;
	}
}

void mvprodc(double complex * vout , double complex * B , double complex * vcolumn , int dim){
	//OBS.: same as mvprodr, but B is a complex
	//computes vout  <- B . vcolumn
	//B is a complex (dim x dim) matrix
	//vcolumn is a complex (dim x 1) column vector
	//vout is a complex (dim x 1) column vector
	int i,j;
	double complex * vtemp = malloc(sizeof(double complex)*dim);

	for(i=0 ; i<dim ; i++){
		vtemp[i] = 0E0;
		for(j=0 ; j<dim ; j++){
			vtemp[i] +=  *getelmc(B,i,j,dim)*vcolumn[j];
		}
	}

	for(i=0 ; i<dim ; i++){
		vout[i] = vtemp[i];
	}
	free(vtemp);
}

void vmprodc(double complex * vout , double complex * vline , double complex * B , int dim){
	//OBS.: same as mvprodr, but B is a complex
	//computes vout  <- vline . B
	//B is a complex (dim x dim) matrix and
	//vline is a complex (1 x dim) line vector
	//vout is a complex (1 x dim) line vector
	int i,j;
	double complex * vtemp = malloc(sizeof(double complex)*dim);

	for(i=0 ; i<dim ; i++){
		vtemp[i] = 0E0;
		for(j=0 ; j<dim ; j++){
			vtemp[i] +=  *getelmc(B,i,j,dim)*vline[j];
		}
	}

	for(i=0 ; i<dim ; i++){
		vout[i] = vtemp[i];
	}
	free(vtemp);
}

double complex det2mc(double complex * v){
	return(v[0]*v[3] - v[1]*v[2]);
}

/////////////////////////////// LINEAR SYSTEMS ALGORITHMS

void bicgmr(double complex * x , double * A , double complex * b , int dim){
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a real (dim x dim) matrix (it doesn't need to be real)

	double tol = 1E-13;
	int icont = 0;
	double complex alpha, beta, auxr;
	double complex * auxv = malloc(sizeof(double complex)*dim);
	double complex * r = malloc(sizeof(double complex)*dim);
	double complex * rbar = malloc(sizeof(double complex)*dim);
	double complex * p = malloc(sizeof(double complex)*dim);
	double complex * pbar = malloc(sizeof(double complex)*dim);

	double * At = malloc(sizeof(double)*dim*dim);
	transpmr(At,A,dim);
	//chute inicial passado em x
	//definições iniciais:

	mvprodr(auxv , A , x , dim);
	cprodvc(auxv , -1.0 , auxv , dim);
	sumvc(r , b , auxv , dim);				//r <- b - A.x


	copyvc(rbar,r,dim);						//rbar <- r
	copyvc(p , r,dim);						//p <- r
	copyvc(pbar , rbar,dim);				//pbar <- rbar

	//printvc(x,dim);
	//printvc(rbar,dim);
	//printf("|r| = %.2E	|rbar| = %.2E \n",normvc(r,dim),normvc(rbar,dim));
	while( (normvc(r,dim) > tol) && (normvc(rbar,dim) > tol) ){			//while( |r| < tol and |rbar| < tol )
		icont++;
		//getchar();getchar();

		mvprodr(auxv , A , p , dim);

		//printc(inprodvc(pbar , auxv,dim));
		//getchar();

		alpha = inprodvc(rbar,r,dim)/inprodvc(pbar , auxv,dim);		//alpha <- (rbar , r)/(pbar , A.p)

		cprodvc(auxv , alpha , p , dim);
		sumvc(x , x , auxv , dim);									//x <- x + alpha . p (SOLUTION STEP)

		auxr = inprodvc(rbar,r,dim);								//auxr <- (rbar,r) to future use in beta

		mvprodr(auxv , A , p , dim);
		cprodvc(auxv , -alpha , auxv , dim);
		sumvc(r , r , auxv , dim);									//r <- r - alpha . A . p

		mvprodr(auxv , At , pbar , dim);
		cprodvc(auxv , -conj(alpha) , auxv , dim);
		sumvc(rbar , rbar , auxv , dim);							//rbar <- rbar - alpha* . At . pbar

		beta = inprodvc(rbar,r,dim)/auxr;							//beta <- (current)(rbar , r)/(last)(rbar , r)

		cprodvc(auxv , beta , p , dim);
		sumvc(p , r , auxv , dim);									//p <- r + beta . p

		cprodvc(auxv , conj(beta) , pbar , dim);
		sumvc(pbar , rbar , auxv , dim);							//pbar <- rbar + beta* . pbar

		//printf("|r| = %.2E	|rbar| = %.2E \n",normvc(r,dim),normvc(rbar,dim));
	}

	//x is the solution to the given tolerance tol

	printf("\nIt took %d iterations to solve the system\n",icont);
	free(r);
	free(rbar);
	free(p);
	free(pbar);
	free(At);
	free(auxv);
}

void bicgmc(double complex * x , double complex * A , double complex * b , int dim){
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a real (dim x dim) matrix (it doesn't need to be real)

	double tol = 1E-13;
	int icont = 0;
	double complex alpha, beta, auxr;
	double complex * auxv 	= malloc(sizeof(double complex)*dim);
	double complex * r 		= malloc(sizeof(double complex)*dim);
	double complex * rbar 	= malloc(sizeof(double complex)*dim);
	double complex * p 		= malloc(sizeof(double complex)*dim);
	double complex * pbar 	= malloc(sizeof(double complex)*dim);

	double complex * Adagger = malloc(sizeof(double complex)*dim*dim);
	hermcmc(Adagger,A,dim);
	//chute inicial passado em x
	//definições iniciais:

	mvprodc(auxv , A , x , dim);
	cprodvc(auxv , -1.0 , auxv , dim);
	sumvc(r , b , auxv , dim);				//r <- b - A.x


	copyvc(rbar,r,dim);						//rbar <- r
	copyvc(p , r,dim);						//p <- r
	copyvc(pbar , rbar,dim);				//pbar <- rbar

	//printvc(x,dim);
	//printvc(rbar,dim);
	//printf("|r| = %.2E	|rbar| = %.2E \n",normvc(r,dim),normvc(rbar,dim));
	while( (normvc(r,dim) > tol) && (normvc(rbar,dim) > tol) ){			//while( |r| < tol and |rbar| < tol )
		icont++;
		//getchar();

		mvprodc(auxv , A , p , dim);
		alpha = inprodvc(rbar,r,dim)/inprodvc(pbar , auxv,dim);		//alpha <- (rbar , r)/(pbar , A.p)

		cprodvc(auxv , alpha , p , dim);
		sumvc(x , x , auxv , dim);									//x <- x + alpha . p (SOLUTION STEP)

		auxr = inprodvc(rbar,r,dim);								//auxr <- (rbar,r) to future use in beta

		mvprodc(auxv , A , p , dim);
		cprodvc(auxv , -alpha , auxv , dim);
		sumvc(r , r , auxv , dim);									//r <- r - alpha . A . p

		mvprodc(auxv , Adagger , pbar , dim);
		cprodvc(auxv , -conj(alpha) , auxv , dim);
		sumvc(rbar , rbar , auxv , dim);							//rbar <- rbar - alpha* . Adagger . pbar

		beta = inprodvc(rbar,r,dim)/auxr;							//beta <- (current)(rbar , r)/(last)(rbar , r)

		cprodvc(auxv , beta , p , dim);
		sumvc(p , r , auxv , dim);									//p <- r + beta . p

		cprodvc(auxv , conj(beta) , pbar , dim);
		sumvc(pbar , rbar , auxv , dim);							//pbar <- rbar + beta* . pbar

		printf("|r| = %.2E	|rbar| = %.2E \n",normvc(r,dim),normvc(rbar,dim));
	}

	//x is the solution to the given tolerance tol

	printf("\nIt took %d iterations to solve the system\n",icont);
	free(r);
	free(rbar);
	free(p);
	free(pbar);
	free(Adagger);
	free(auxv);
}

void cgmr(double * x , double * A , double * b , int dim){
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are real (dim x 1) column vectors
	//A is a real symmetric (dim x dim) matrix

	double tol = 1E-13;
	int icont = 0;
	double alpha, beta, auxr;
	double * auxv = malloc(sizeof(double)*dim);
	double * r = malloc(sizeof(double)*dim);
	double * p = malloc(sizeof(double)*dim);

	//chute inicial passado em x
	//definições iniciais:

	mvrprod(auxv , A , x , dim);
	cprodvr(auxv , -1.0 , auxv , dim);
	sumvr(r , b , auxv , dim);				//r <- b - A.x

	copyvr(p , r,dim);						//p <- r

	printf("|r| = %.2E\n",normvr(r,dim));
	while( (normvr(r,dim) > tol) ){			//while( |r| < tol )
		icont++;

		mvrprod(auxv , A , p , dim);

		auxr = inprodvr(r,r,dim);									//auxr <- (rbar,r) to future use in beta

		alpha = auxr/inprodvr(p , auxv,dim);						//alpha <- (r , r)/(p , A.p)


		printc(alpha);
		getchar();

		cprodvr(auxv , alpha , p , dim);
		sumvr(x , x , auxv , dim);									//x <- x + alpha . p (SOLUTION STEP)

		mvrprod(auxv , A , p , dim);
		cprodvr(auxv , -alpha , auxv , dim);
		sumvr(r , r , auxv , dim);									//r <- r - alpha . A . p

		beta = inprodvr(r,r,dim)/auxr;							//beta <- (current)(r , r)/(last)(r , r)

		cprodvr(auxv , beta , p , dim);
		sumvr(p , r , auxv , dim);									//p <- r + beta . p

		printf("|r| = %.2E \n",normvr(r,dim));
	}

	//x is the solution to the given tolerance tol

	printf("\nIt took %d iterations to solve the system\n",icont);
	free(r);
	free(p);
	free(auxv);
}

void cgmc(double complex * x , double complex * A , double complex * b , int dim){
	//algorithm based on Numerical Recipes book
	//linear system A . x = b
	//x and b are complex (dim x 1) column vectors
	//A is a real (dim x dim) matrix (it doesn't need to be real)

	double tol = 1E-13;
	int icont = 0;
	double complex alpha, beta, auxr;
	double complex * auxv = malloc(sizeof(double complex)*dim);
	double complex * r = malloc(sizeof(double complex)*dim);
	double complex * p = malloc(sizeof(double complex)*dim);

	//chute inicial passado em x
	//definições iniciais:

	mvprodc(auxv , A , x , dim);
	cprodvc(auxv , -1.0 , auxv , dim);
	sumvc(r , b , auxv , dim);				//r <- b - A.x

	copyvc(p , r,dim);						//p <- r

	//printvc(x,dim);
	//printvc(rbar,dim);
	//printf("|r| = %.2E\n",normvc(r,dim));
	while( (normvc(r,dim) > tol) ){			//while( |r| < tol )
		icont++;

		mvprodc(auxv , A , p , dim);

		auxr = inprodvc(r,r,dim);									//auxr <- (r,r) to future use in beta
		alpha = auxr/inprodvc(p , auxv,dim);						//alpha <- (r , r)/(p , A.p)

		cprodvc(auxv , alpha , p , dim);
		sumvc(x , x , auxv , dim);									//x <- x + alpha . p (SOLUTION STEP)

		mvprodc(auxv , A , p , dim);
		cprodvc(auxv , -alpha , auxv , dim);
		sumvc(r , r , auxv , dim);									//r <- r - alpha . A . p

		beta = inprodvc(r,r,dim)/auxr;							//beta <- (current)(r , r)/(last)(r , r)

		cprodvc(auxv , beta , p , dim);
		sumvc(p , r , auxv , dim);									//p <- r + beta . p

		//printf("|r| = %.2E \n",normvc(r,dim));
	}

	//x is the solution to the given tolerance tol

	printf("\nIt took %d iterations to solve the system\n",icont);
	free(r);
	free(p);
	free(auxv);
}

///////////////////////////////

void setzeroqvector(int * v){
		v[0]=0;
		v[1]=0;
		v[2]=0;
		v[3]=0;
	}

int aseps(int i , int j , int k){
	i++;
	j++;
	k++;
	if(i==j || j==k || i==k){
		return(0);
	}
	else if( (i==1 && j==2 && k==3) || (i==3 && j==1 && k==2) || (i==2 && j==3 && k==1)){
		return(1);
	}
	else{
		return(-1);
	}
}

//////////////////////////////////////NOT INCORPORATED YET

void setonevc(double complex * vc,int dim){
	//set all elements of complex vector to unity
	int i;
	for(i=0;i<dim;i++){
		vc[i] = 1e0;
	}
}

void setrandomvc(double complex * vc,int dim){
	//set all elements of complex vector to a random number
	int i;
	for(i=0;i<dim;i++){
		vc[i] = ran0(global_seed)/sqrt(dim);
	}
}

void setonevr(double * vr,int dim){
	//set all elements of complex vector to unity
	int i;
	for(i=0;i<dim;i++){
		vr[i] = 1e0;
	}
}

void setplanewave(double complex * target, double * p , int atarget){
	//set target vector to plane wave to exp(-2*pi*i*(p.x)) within
	//the subspace of a==atarget and zero to the remaining
	int a,t,x,y,z;
	for(a=0;a<3;a++){
		for(t=0;t<Nt;t++){
			for(x=0;x<N;x++){
				for(y=0;y<N;y++){
					for(z=0;z<N;z++){
						if(a==atarget){
							*getelvc(target,a,t,x,y,z) = cexp( -2.0*M_PI*I*(p[0]*t + p[1]*x + p[2]*y + p[3]*z) );
						}
						else{
							*getelvc(target,a,t,x,y,z) = 0E0;
						}
					}
				}
			}
		}
	}
}

void setplanewaveallcolors(double complex * target, double * p ){
	//set target vector to plane wave to exp(-2*pi*i*(p.x)) within
	//the subspace of a==atarget and zero to the remaining
	int a,t,x,y,z;
	for(a=0;a<3;a++){
		for(t=0;t<Nt;t++){
			for(x=0;x<N;x++){
				for(y=0;y<N;y++){
					for(z=0;z<N;z++){
							*getelvc(target,a,t,x,y,z) = cexp( -2.0*M_PI*I*(p[0]*t + p[1]*x + p[2]*y + p[3]*z) );
					}
				}
			}
		}
	}
}

void rsetplanewave(double * target, double * p , int atarget){
	//set target vector to plane wave to exp(-2*pi*i*(p.x)) within
	//the subspace of a==atarget and zero to the remaining
	int a,t,x,y,z;
	for(a=0;a<3;a++){
		for(t=0;t<Nt;t++){
			for(x=0;x<N;x++){
				for(y=0;y<N;y++){
					for(z=0;z<N;z++){
						if(a==atarget){
							*getelvr(target,a,t,x,y,z) = creal(cexp( -2.0*M_PI*I*(p[0]*t + p[1]*x + p[2]*y + p[3]*z) ));
						}
						else{
							*getelvr(target,a,t,x,y,z) = 0E0;
						}
					}
				}
			}
		}
	}
}

void rsetplanewaveallcolors(double * target, double * p){
	//set target vector to plane wave to exp(-2*pi*i*(p.x)) within
	//the subspace of a==atarget and zero to the remaining
	int a,t,x,y,z;
	for(a=0;a<3;a++){
		for(t=0;t<Nt;t++){
			for(x=0;x<N;x++){
				for(y=0;y<N;y++){
					for(z=0;z<N;z++){
							*getelvr(target,a,t,x,y,z) = creal(cexp( -2.0*M_PI*I*(p[0]*t + p[1]*x + p[2]*y + p[3]*z) ));
					}
				}
			}
		}
	}
}

void setrandomvr(double * vr,int dim){
	//set all elements of real vector to a random number
	int i;
	for(i=0;i<dim;i++){
		vr[i] = ran0(global_seed)/sqrt(dim);
	}
}

void invMv(double * vout, double * vin){
	double complex * temp = malloc(sizeof(double complex)*4);
	double det = detv(vin);
	temp[0] = getelm(vin,1,1);
	temp[1] = -getelm(vin,0,1);
	temp[2] = -getelm(vin,1,0);
	temp[3] = getelm(vin,0,0);

	cprodvc(temp,1/det,temp,2);

	converttov(vout,temp);
}

void invMc(double complex * vout, double complex * vin){
	double complex * temp = malloc(sizeof(double complex)*4);
	double complex det = det2mc(vin);
	temp[0] = vin[3];
	temp[1] = -vin[1];
	temp[2] = -vin[2];
	temp[3] = vin[0];

	cmprodc(temp,1e0/det,temp,2);
	copymc(vout,temp,2);
}

int eps4(int a, int b , int c , int d){
	int bubbleCount4(int v1,int v2,int v3,int v4){
		int pcount=0, i,j, aux;
		int v[4] = {a,b,c,d};
		int size = 4;
		for(i=1;i<size;i++){
			for(j=0;j<size-i;j++){
				if(v[j]>v[j+1]){
					aux = v[j];
					v[j] = v[j+1];
					v[j+1] = aux;
					pcount++;
				}
			}
		}
		return(pcount);
	}
	if( (a==b)||(a==c)||(a==d)||(b==c)||(b==d)||(c==d) ){
		return(0);
	}
	else if( (bubbleCount4(a,b,c,d) % 2) == 0 ){			//is even permutation
		return(1);
	}
	else{
		return(-1);
	}
}

int eps4tilde(int a, int b, int c, int d){
	return( signal(a*b*c*d)*eps4(abs(a),abs(b),abs(c),abs(d)) );
}
