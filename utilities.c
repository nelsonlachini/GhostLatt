#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "utilities.h"
#include "global.h"
#include "algebra.h"

////////////////////////////// ran0 implementation from numerical recipes
static const int IA = 16807;
static const long IM = 2147483647;
static const double AM = (1.0/IM);
static const int IQ = 127773;
static const int IR = 2836;
static const long MASK = 123459876;

float ran0(long *idum)
{	
	long k;
	float ans;
	*idum ^= MASK;
	k = (*idum)/IQ;
	*idum = IA*(*idum - k*IQ) - IR*k;
	if (*idum < 0) *idum += IM;
	ans = AM*(*idum);
	*idum ^= MASK;
	return ans;
}

////////////////////////////// SU2 MATRIX/VECTORS METHODS

double * getelv(double * vector , int i){
	return vector + i;
}

double complex getelm(double * vector , int i, int j){
	if(i==0){
		if(j==0){
			return *getelv(vector , 0) + I*(*getelv(vector , 3));
		}
		else{
			return *getelv(vector , 2) + I*(*getelv(vector , 1));
		}
	}
	else{
		if(j==0){
			return -*getelv(vector , 2) + I*(*getelv(vector , 1));
		}
		else{
			return *getelv(vector , 0) - I*(*getelv(vector , 3));
		}
	}
}

void printv(double * vector){
	printf("(%lf , %lf , %lf , %lf) \n", vector[0],vector[1],vector[2],vector[3]);
}

void printc(double complex z){
	//function to print a complex doubles
	printf(" %e + %e i ", creal(z) , cimag(z) );
}

void printr(double r){
	//function to print a doubles
	printf(" %e ", r );
}

void fprintc( FILE * file , double complex z){
	fprintf(file, " %lf + %lf i", creal(z) , cimag(z) );
}

void printm(double * A){
	//function to print a square matrix (2x2) of complex doubles
	int i,j;
	for(i=0 ; i<2 ; i++){
		for(j=0 ; j<2 ; j++){
			printc( getelm(A,i,j) );
		}
		printf("\n");
	}
	printf("\n");
}

void fprintm(FILE * file , double * A){
	int i,j;
	for(i=0 ; i<2 ; i++){
		for(j=0 ; j<2 ; j++){
			fprintc(file,  getelm(A,i,j) );
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

double * getU(double * lattice, int t , int x , int y, int z , int mi){
	return ( lattice + t*spatialV*16 + x*N2*16 + y*N*16 + z*16 + mi*4);
}

double * getUdim(double * lattice, int t , int x , int y, int z , int mi, int N){
	return ( lattice + t*spatialV*16 + x*N2*16 + y*N*16 + z*16 + mi*4);
}

void printLattice(double * lattice){
	int t,x,y,z,mi;
	printf("------------------");

	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						printf("\n");
						printm( getU(lattice,t,x,y,z,mi) );
					}
				}
			}
		}
	}
	printf("------------------");
}

void copyv(double * out , double * in){
	//copy the values of the elements of in matrix to the out matrix
	int i;
	for(i=0 ; i<4 ; i++){
		*( getelv(out,i) ) = *( getelv(in,i) );
	}
}

double * randSU2v(double * vector ,long * seed , double parameter){
	int j;
	float eps = parameter-0.01;				//displacement from identity
	double r[3] , x[4];

	x[0] = signal( ran0(seed) - 0.5 )*sqrt(1-eps*eps);

	for(j=0 ; j < 3 ; j++){
		r[j] = ran0(seed) - 0.5;
	}

	for(j=1 ; j < 4 ; j++){
		x[j] = eps*r[j-1]/sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
	}

	for(j=0 ; j < 4 ; j++){
		*getelv(vector , j) = x[j];
	}
}


////////////////////////////// GENERAL MATRIX FUNCTIONS WITH DOUBLE VARIABLES

double * getelmr(double * matrix , int i , int j , int dim){
	return(matrix + i*dim + j);
}

void printmr(double * A , int dim){
	//function to print a square matrix (dim x dim) of doubles
	int i,j;
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			printf(" %E ", *getelmr(A,i,j,dim) );
		}
		printf("\n");
	}
	printf("\n");
}

void copymr(double * Mout , double * Min , int dim){
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmr(Mout,i,j,dim) = *getelmr(Min,i,j,dim);
		}
	}
}

void copyvr(double * vout, double * vin , int dim){
	int i;
	for(i=0;i<dim;i++){
		vout[i] = vin[i];
	}
}

double * getelvr(double * vr , int a , int t , int x , int y , int z){
	return(vr + N*(N*(N*(N*a+t) +x) + y) + z);
}

void printvr(double * V , int dim){
	//function to print a (dim x 1) vector of complex variables
	int i;
	printf("\n-----------------");
	for(i=0;i<dim;i++){
		printf("\n");
		printr( V[i] );
	}
	printf("\n-----------------");
}

////////////////////////////// GENERAL MATRIX FUNCTIONS WITH DOUBLE COMPLEX VARIABLES

double complex * getelmc(double complex * M , int i , int j , int dim){
	return(M + i*dim + j);
}

void printvc(double complex * V , int dim){
	//function to print a (dim x 1) vector of complex variables
	int i;
	printf("\n-----------------");
	for(i=0;i<dim;i++){
		printf("\n");
		printc( V[i] );
	}
	printf("\n-----------------");
}

void copyvc(double complex * Vout , double complex * Vin, int dim){
	int i;
	for(i=0;i<dim;i++){
		Vout[i] = Vin[i];
	}
}

void copymc(double complex * Mout , double complex * Min , int dim){
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			*getelmc(Mout,i,j,dim) = *getelmc(Min,i,j,dim);
		}
	}
}

void printmc(double complex * A , int dim){
	//function to print a square matrix (dim x dim) of complex doubles
	int i,j;
	for(i=0 ; i<dim ; i++){
		for(j=0 ; j<dim ; j++){
			printc( *getelmc(A,i,j,dim) );
		}
		printf("\n");
	}
	printf("\n");
}

double complex * getelvc(double complex * vc , int a , int t , int x , int y , int z){
	return(vc + a*totalV + t*spatialV + x*N2 + y*N + z);
}

void converttov(double * vout, double complex * vin){
	vout[0] = creal(*getelmc(vin,0,0,2));
	vout[1] = cimag(*getelmc(vin,0,1,2));
	vout[2] = creal(*getelmc(vin,0,1,2));
	vout[3] = cimag(*getelmc(vin,0,0,2));
}

////////////////////////////// OTHER STUFF

void setquadv(int * amu, int mu){
	int i;
	for(i=0;i<4;i++){
		if(abs(mu)==i){
			amu[i] = 1*signal(mu);
		}
		else{
			amu[i] = 0;
		}
	}
}

void initl(double * lattice , float parameter){
	int t,x,y,z,mi,a,b;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						if(parameter == 0){
							setidv( getU(lattice, t,x,y,z,mi) );
						}
						else{
							randSU2v(getU(lattice, t,x,y,z,mi) , global_seed , parameter);
						}
					}
				}
			}
		}
	}
}

void reunitl(double * lattice){
	int t,x,y,z,mi,i;
	double * aux_link;
	double norm;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						aux_link = getU(lattice, t,x,y,z,mi);
						norm = sqrt(detv(aux_link));
						//printc(detv(aux_link));
						for(i=0;i<4;i++){
							*getelv(aux_link , i) /= norm;
						}
						//printr(detv(aux_link));
						//printf("\n");
					}
				}
			}
		}
	}
}

void reunitv(double * vector){
	int i;
	double norm;
	norm = sqrt(detv(vector));
	for(i=0;i<4;i++){
		*getelv(vector , i) /= norm;
	}
}

double modulus(double a, double b){
		int result = (int)( a / b );
		return a - (double)( result ) * b;
	}

double getStep( double coord , double dcoord ){
	//return the lattice coordinate given the increment(useful when there is a chance to get a negative coordinate)
	double aux = modulus(coord+dcoord,N);
	if(aux<0)
		return N+aux;
	return aux;
}

double getStepT( double coord , double dcoord ){
	//return the lattice coordinate given the increment(useful when there is a chance to get a negative coordinate)
	double aux = modulus(coord+dcoord,Nt);
	if(aux<0)
		return Nt+aux;
	return aux;
}

void copyl(double * Lout , double * Lin){
	int t,x,y,z,mi,a,b;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						copyv( getU(Lout, t,x,y,z,mi) , getU(Lin,t,x,y,z,mi) );
					}
				}
			}
		}
	}
}

double * getg(double * g, int t , int x , int y, int z){
	return ( g + t*spatialV*4 + x*N2*4 + y*N*4 + z*4);
}

void printvcbycolor(double complex * v, int a){
	int t,x,y,z;
	for(a=0;a<3;a++){
		for(t=0;t<Nt;t++){
			for(x=0;x<N;x++){
				for(y=0;y<N;y++){
					for(z=0;z<N;z++){
						printf("%d %d %d %d %d",a,t,x,y,z);
						printc(*getelvc(v,a,t,x,y,z));
						printf("\n");
					}
				}
			}
		}
	}
}

void tseed(){
	//seed based in the system time
	global_seed = malloc(sizeof(long));
	*global_seed = time(0);
}

void progress_panel(int i, int total){
		printf("\rIn progress... %.1f%%",(100.0*i)/total);
		fflush(stdout);
	}

void save_lattice(double * lattice, char * file_name){
	//always use double quotes to pass file_name
	FILE * fl;
	fl = fopen( file_name , "wb");
	if( fwrite(lattice, sizeof(double) , dimLattice , fl) ) {
		printf("\n Lattice saved");
	}
	else{
		perror("Error saving lattice");
		exit(EXIT_FAILURE);
	}
	fclose(fl);
}

void save_vector(double complex * vector, char * file_name){
	//always use double quotes to pass file_name
	FILE * fl;
	fl = fopen( file_name , "wb");
	if( fwrite(vector, sizeof(double complex) , colorV , fl) ) {
		//printf("\n Vector saved");
	}
	else{
		perror("Error saving vector");
		exit(EXIT_FAILURE);
	}
	fclose(fl);
}

void save_rvector(double * vector, char * file_name){
	//always use double quotes to pass file_name
	FILE * fl;
	fl = fopen( file_name , "wb");
	if( fwrite(vector, sizeof(double) , colorV , fl) ) {
		//printf("\n Vector saved");
	}
	else{
		perror("Error saving vector");
		exit(EXIT_FAILURE);
	}
	fclose(fl);
}

void getName(char * out, char * header , double num){
	//generates string name : header+num
	//always use double quotes to pass header
	char aux[15];
	strcpy(out, "");
	strcpy(out , header);
	if( num - (int)num == 0){
		sprintf(aux , "%d" , (int)num);
	}
	else{
		sprintf(aux , "%.2lf" , num);
	}
	strcat(out , aux);

}

void load_lattice(double * lattice , char * file_name){
	FILE * fl;
	fl = fopen( file_name, "rb" );
	if( fread( lattice, sizeof(double) , dimLattice , fl) ){
		printf("\n Lattice loaded");
	}
	else{
		perror("Error loading lattice");
		exit(EXIT_FAILURE);
	fclose(fl);
	}
}

void load_rvector(double * vector , char * file_name){
	FILE * fl;
	fl = fopen( file_name, "rb" );
	if( fread( vector, sizeof(double) , colorV , fl) ){
		printf("\n Vector loaded");
	}
	else{
		perror("Error loading vector");
		exit(EXIT_FAILURE);
	}
	fclose(fl);
}

int compareLink(double * link1, double * link2){
	int i;
	double tol = 1e-8;
	for(i=0;i<4;i++){
		if( fabs(link1[i] - link2[i]) > tol ){
			//printc(fabs(link1[i]-link2[i]));
			return 0;
		}
	}
	return 1;
}

int compareLattice(double * lattice1, double * lattice2){
	int t,x,y,z,mi,a,b;
	double complex * aux = malloc(sizeof(double complex)*4);
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						if( !compareLink( getU(lattice1, t,x,y,z,mi) , getU(lattice2,t,x,y,z,mi) )  ){
							printf("\n ----> Different lattices!");
							return 0;
						}
					}
				}
			}
		}
	}
	printf("\n ----> Same lattices!");
	return 1;
}

int comparevc(double complex * v1, double complex * v2, int dim){
	int i;
	double tol = 1e-13;
	for(i=0;i<dim;i++){
		if( cabs(v1[i] - v2[i]) > tol){
			printf("\n ----> Different vectors! <----");
			return(0);
		}
	}
	printf("\n ----> Same vectors! <----");
	return(1);
}

int comparevr(double * v1, double * v2, int dim){
	int i;
	double tol = 1e-13;
	for(i=0;i<dim;i++){
		if( cabs(v1[i] - v2[i]) > tol){
			printf("\n ----> Different vectors! <----");
			return(0);
		}
	}
	printf("\n ----> Same vectors! <----");
	return(1);
}

void reescaleGaugeField(double * out, double * in, double scale){
	int t,x,y,z,mi,a;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						getU(out,t,x,y,z,mi)[0] = getU(in,t,x,y,z,mi)[0];
						for(a=1;a<4;a++)
							getU(out,t,x,y,z,mi)[a] = getU(in,t,x,y,z,mi)[a]*scale;
					}
				}
			}
		}
	}

}

void reescaleLinks(double * out, double * in, double scale){
	int t,x,y,z,mi,a;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4 ; mi++){
						getU(out,t,x,y,z,mi)[0] = getU(in,t,x,y,z,mi)[0]*scale;
					}
				}
			}
		}
	}

}

double compareVectorToPW(double * eigenvector, double pPW){
		double totalsum = 0e0;
		double complex fouriersum = 0e0;
		double complex pPW2ipi = -2*M_PI*I*pPW;
		int x[4];
		int b,mu;

		for(b=0;b<3;b++){
			for(mu=0;mu<4;mu++){
				fouriersum=0;
				for(x[0]=0;x[0]<Nt;x[0]++){
					for(x[1]=0;x[1]<N;x[1]++){
						for(x[2]=0;x[2]<N;x[2]++){
							for(x[3]=0;x[3]<N;x[3]++){
								fouriersum += *getelvr(eigenvector,b,x[0],x[1],x[2],x[3])*cexp(pPW2ipi*x[mu]) ;
							}
						}
					}
				}
				totalsum += pow( cabs(fouriersum),2 );
			}
		}

		return(totalsum/(3.0*4.0*totalV));
}

double zeroMomentumTransform(double * vector){
	double fouriersum = 0e0;
	int x[4];
	int b;

	for(b=0;b<3;b++){
		fouriersum=0e0;
		for(x[0]=0;x[0]<Nt;x[0]++){
			for(x[1]=0;x[1]<N;x[1]++){
				for(x[2]=0;x[2]<N;x[2]++){
					for(x[3]=0;x[3]<N;x[3]++){
						fouriersum += *getelvr(vector,b,x[0],x[1],x[2],x[3]);
					}
				}
			}
		}
	}

	return(fouriersum/(3e0*totalV));
}

int isLatticeUnitary(double * lattice){
	int t,x,y,z,mi,i;
	double * aux_link;
	double norm;
	for(t=0 ; t<Nt ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mi=0 ; mi<4; mi++){
						if( !(fabs(detv(getU(lattice,t,x,y,z,mi))-1e0) < 1e-4) ){
							printf("\n Lattice is not unitarized.");
							return(0);
						}
					}
				}
			}
		}
	}
	printf("\n Lattice is unitarized.");
	return(1);
}

int isLinkUnitary(double * link){
	if( fabs(detv(link)-1e0) < 1e-4 ){
		printf("\n Link is unitarized.");
		return(1);
	}
	printf("\n Link is not unitarized.");
	return(0);
}

int isLatticeNaN(double * lattice){
	int i;
	for(i=0;i<dimLattice;i++){
		if( isnan(lattice[dimLattice]) ){
			//printf("\n WARNING! LATTICE IS NAN!\n");
			return(1);
		}
	}
	//printf("\n WARNING! LATTICE IS NOT NAN!\n");
	return(0);
}

double max(double * array, int dim){
	double maximum = array[0];
	int i;
	for(i=1;i<dim;i++)
		if(array[i]>maximum)
			maximum=array[i];
	return(maximum);
}

double min(double * array, int dim){
	double minimum = array[0];
	int i;
	for(i=1;i<dim;i++)
		if(fabs(array[i])<fabs(minimum))
			minimum=array[i];
	return(minimum);
}

double sortDouble(double i , double j){
    //sort int in the interval [i,j)
    return((double)(ran0(global_seed)*(j-i) + i));
}

void printpos(double t, double x, double y , double z, double mi){
	printf("\n %lf | %lf | %lf | %lf | %lf \n",t,x,y,z,mi);
}