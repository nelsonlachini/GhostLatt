#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "utilities.h"
#include "algebra.h"
#include "global.h"

////////////////////////////// ran0 implementation from numerical recipes
float ran0(long *idum)
{	
	static const int IA = 16807;
	static const long IM = 2147483647;
	static const double AM = (1.0/IM);
	static const int IQ = 127773;
	static const int IR = 2836;
	static const long MASK = 123459876;
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

double * getLink(LatticeLinkSU2 * lattice, int t , int x , int y, int z , int mi){
	return ( lattice->U + t*lattice->spatialV*16 + x*lattice->N2*16 + y*lattice->N*16 + z*16 + mi*4);
}

// double * getLinkdim(LatticeLinkSU2 * lattice, int t , int x , int y, int z , int mi, int N){
// 	return ( lattice->U + t*lattice->spatialV*16 + x*lattice->N2*16 + y*lattice->N*16 + z*16 + mi*4);
// }

void printLattice(LatticeLinkSU2 * lattice){
	int t,x,y,z,mi;
	printf("------------------");
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4; mi++){
			printf("\n");
			printm( getLink(lattice,t,x,y,z,mi) );
		}
	}}}}
	printf("------------------");
}

void copyv(double * out , double * in){
	//copy the values of the elements of in matrix to the out matrix
	int i;
	for(i=0 ; i<4 ; i++){
		*( getelv(out,i) ) = *( getelv(in,i) );
	}
}

void randSU2v(double * vector ,long * seed , double parameter){
	int j;
	float eps = parameter-0.01;				//displacement from identity
	double r[3] , x[4];

	x[0] = getSignal( ran0(seed) - 0.5 )*sqrt(1-eps*eps);

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

void copyvr(double * vout, double * vin, int dim){
	int i;
	for(i=0 ; i<dim ; i++)
		vout[i] = vin[i];
}

void copyLatticeColorVectorReal(LatticeColorVectorReal * vout, LatticeColorVectorReal * vin){
	int i;
	//check dimensions
	if((vout->N != vin->N) || (vout->Nt != vin->Nt))
		printf("\nWARNING! Vectors have different dimensions!\n");
	else
		for(i=0 ; i<vout->dimVector ; i++)
			vout->vec[i] = vin->vec[i];
}

double * getelvr(double * vr , int a , int t , int x , int y , int z){
	return(vr + N*(N*(N*(N*a+t) +x) + y) + z);
}

double * getLatticeColorVectorReal(LatticeColorVectorReal * vr , int a , int t , int x , int y , int z){
	return(vr->vec + vr->N*(vr->N*(vr->N*(vr->Nt*a+t) +x) + y) + z );
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

void copyLatticeColorVectorComplex(LatticeColorVectorComplex * vout , LatticeColorVectorComplex * vin){
	int i;
	if((vout->N != vin->N) || (vout->Nt != vin->Nt))
		printf("\nWARNING! Vectors have different dimensions!\n");
	else
		for(i=0;i<vout->dimVector;i++)
			vout->vec[i] = vin->vec[i];
		
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

double complex * getLatticeColorVectorComplex(LatticeColorVectorComplex * vc , int a , int t , int x , int y , int z){
	return(vc->vec + vc->N*(vc->N*(vc->N*(vc->Nt*a+t) +x) + y) + z );
}

void converttov(double * vout, double complex * vin){
	vout[0] = creal(*getelmc(vin,0,0,2));
	vout[1] = cimag(*getelmc(vin,0,1,2));
	vout[2] = creal(*getelmc(vin,0,1,2));
	vout[3] = cimag(*getelmc(vin,0,0,2));
}

////////////////////////////// OTHER STUFF

void setUnitVector(int * amu, int mu){
	int i;
	for(i=0;i<4;i++){
		if(abs(mu)==i){
			amu[i] = 1*getSignal(mu);
		}
		else{
			amu[i] = 0;
		}
	}
}

void initl(LatticeLinkSU2 * lattice , float parameter){
	int t,x,y,z,mi,a,b;
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4; mi++){
			if(parameter == 0){
				setidv( getLink(lattice, t,x,y,z,mi) );
			}
			else{
				randSU2v(getLink(lattice, t,x,y,z,mi) , global_seed , parameter);
			}
		}
	}}}}
}

void reunitLatticeLinkSU2(LatticeLinkSU2 * lattice){
	int t,x,y,z,mi,i;
	double * aux_link;
	double norm;
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4; mi++){
			aux_link = getLink(lattice, t,x,y,z,mi);
			norm = sqrt(detv(aux_link));
			//printc(detv(aux_link));
			for(i=0;i<4;i++){
				*getelv(aux_link , i) /= norm;
			}
			//printr(detv(aux_link));
			//printf("\n");
		}
	}}}}
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

void copyLatticeLinkSU2(LatticeLinkSU2 * Lout , LatticeLinkSU2 * Lin){
	int t,x,y,z,mi,a,b;
	for(t=0 ; t<Lout->Nt ; t++){
	for(x=0 ; x<Lout->N ; x++){
	for(y=0 ; y<Lout->N ; y++){
	for(z=0 ; z<Lout->N ; z++){
		for(mi=0 ; mi<4; mi++){
			copyv( getLink(Lout, t,x,y,z,mi) , getLink(Lin,t,x,y,z,mi) );
		}
	}}}}
}

double * getg(LatticeGaugeSU2 * gauge_transf, int t , int x , int y, int z){
	return ( gauge_transf->g + t*gauge_transf->spatialV*4 + x*gauge_transf->N2*4 + y*gauge_transf->N*4 + z*4);
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
		}}}}
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

void saveLatticeLinkSU2(LatticeLinkSU2 * lattice, char * file_name){
	//always use double quotes to pass file_name
	FILE * fl;
	fl = fopen( file_name , "wb");
	if( fwrite(lattice->U, sizeof(double) , lattice->dimLattice , fl) ) {
		printf("\n Lattice saved");
	}
	else{
		perror("Error saving lattice");
		exit(EXIT_FAILURE);
	}
	fclose(fl);
}

void saveLatticeColorVectorComplex(LatticeColorVectorComplex * vector, char * file_name){
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

void saveLatticeColorVectorReal(LatticeColorVectorReal * vector, char * file_name){
	//always use double quotes to pass file_name
	FILE * fl;
	fl = fopen( file_name , "wb");
	if( fwrite(vector->vec, sizeof(double) , vector->dimVector , fl) ) {
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

void loadLatticeLinkSU2(LatticeLinkSU2 * lattice , char * file_name){
	FILE * fl;
	fl = fopen( file_name, "rb" );
	if( fread( lattice->U, sizeof(double) , lattice->dimLattice , fl) ){
		printf("\n Lattice loaded");
	}
	else{
		perror("Error loading lattice");
		exit(EXIT_FAILURE);
	fclose(fl);
	}
}

void loadLatticeColorVectorReal(LatticeColorVectorReal * vector , char * file_name){
	FILE * fl;
	fl = fopen(file_name, "rb");
	if( fread(vector->vec, sizeof(double), vector->dimVector, fl) ){
		printf("\n Vector loaded");
	}
	else{
		perror("Error loading vector");
		exit(EXIT_FAILURE);
	}
	fclose(fl);
}

void loadLatticeColorVectorComplex(LatticeColorVectorComplex * vector , char * file_name){
	FILE * fl;
	fl = fopen(file_name, "rb");
	if( fread(vector->vec, sizeof(double complex), vector->dimVector, fl) ){
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

int compareLattice(LatticeLinkSU2 * lattice1, LatticeLinkSU2 * lattice2){
	int t,x,y,z,mi,a,b;
	double complex * aux = malloc(sizeof(double complex)*4);
		for(t=0 ; t<lattice1->Nt ; t++){
		for(x=0 ; x<lattice1->N ; x++){
		for(y=0 ; y<lattice1->N ; y++){
		for(z=0 ; z<lattice1->N ; z++){
			for(mi=0 ; mi<4; mi++){
				if( !compareLink( getLink(lattice1, t,x,y,z,mi) , getLink(lattice2,t,x,y,z,mi) )  ){
					printf("\n ----> Different lattices!");
					return 0;
				}
			}
		}}}}
	printf("\n ----> Same lattices!");
	return 1;
}

int compareLatticeColorVectorComplex(LatticeColorVectorComplex * v1, LatticeColorVectorComplex * v2){
	int i;
	double tol = 1e-13;
	if(v1->dimVector != v2->dimVector){
		printf("\n Vectors have different dimensions!");
		return(0);
	}
	for(i=0;i<v1->dimVector;i++){
		if( cabs(v1->vec[i] - v2->vec[i]) > tol){
			printf("\n ----> Different vectors! <----");
			return(0);
		}
	}
	printf("\n ----> Same vectors! <----");
	return(1);
}

int compareLatticeColorVectorReal(LatticeColorVectorReal * v1, LatticeColorVectorReal * v2){
	int i;
	double tol = 1e-13;
	for(i=0;i<v1->dimVector;i++){
		if( fabs(v1->vec[i] - v2->vec[i]) > tol){
			printf("\n ----> Different vectors! <----");
			return(0);
		}
	}
	printf("\n ----> Same vectors! <----");
	return(1);
}

void reescaleGaugeField(LatticeLinkSU2 * out, LatticeLinkSU2 * in, double scale){
	int t,x,y,z,mi,a;
	for(t=0 ; t<out->Nt ; t++){
	for(x=0 ; x<out->N ; x++){
	for(y=0 ; y<out->N ; y++){
	for(z=0 ; z<out->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			getLink(out,t,x,y,z,mi)[0] = getLink(in,t,x,y,z,mi)[0];
			for(a=1;a<=3;a++)
				getLink(out,t,x,y,z,mi)[a] = getLink(in,t,x,y,z,mi)[a]*scale;
		}
	}}}}
}

void reescaleLinks(LatticeLinkSU2 * out , LatticeLinkSU2 * in, double scale){
	int t,x,y,z,mi,a;
	for(t=0 ; t<out->Nt ; t++){
	for(x=0 ; x<out->N ; x++){
	for(y=0 ; y<out->N ; y++){
	for(z=0 ; z<out->N ; z++){
		for(mi=0 ; mi<4 ; mi++){
			getLink(out,t,x,y,z,mi)[0] = getLink(in,t,x,y,z,mi)[0]*scale;
		}
	}}}}
}

double compareVectorToPW(LatticeColorVectorReal * eigenvector, double pPW){
		double totalsum = 0e0;
		double complex fouriersum = 0e0;
		double complex pPW2ipi = -2*M_PI*I*pPW;
		int x[4];
		int b,mu;

		for(b=0;b<3;b++){
			for(mu=0;mu<4;mu++){
				fouriersum=0;
				for(x[0]=0;x[0]<eigenvector->Nt;x[0]++){
				for(x[1]=0;x[1]<eigenvector->N;x[1]++){
				for(x[2]=0;x[2]<eigenvector->N;x[2]++){
				for(x[3]=0;x[3]<eigenvector->N;x[3]++){
					fouriersum += *getLatticeColorVectorReal(eigenvector,b,x[0],x[1],x[2],x[3])*cexp(pPW2ipi*x[mu]) ;
				}}}}
				totalsum += pow( cabs(fouriersum),2 );
			}
		}

		return(totalsum/(4.0*eigenvector->dimVector));
}

/*
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
		}}}}
	}
	return(fouriersum/(3e0*totalV));
}*/

int isLatticeUnitary(LatticeLinkSU2 * lattice){
	int t,x,y,z,mi,i;
	double * aux_link;
	double norm;
	for(t=0 ; t<lattice->Nt ; t++){
	for(x=0 ; x<lattice->N ; x++){
	for(y=0 ; y<lattice->N ; y++){
	for(z=0 ; z<lattice->N ; z++){
		for(mi=0 ; mi<4; mi++){
			if( !(fabs(detv(getLink(lattice,t,x,y,z,mi))-1e0) < 1e-4) ){
				printf("\n Lattice is not unitarized.");
				return(0);
			}
		}
	}}}}
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

int isLatticeNaN(LatticeLinkSU2 * lattice){
	int i;
	for(i=0;i<lattice->dimLattice;i++){
		if( isnan(lattice->U[i]) ){
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

LatticeLinkSU2* newLatticeLinkSU2(int _N, int _Nt) {    // emulates a constructor for lattice object 
	LatticeLinkSU2* lattice_out = malloc(sizeof(LatticeLinkSU2));
	lattice_out->N = _N;
	lattice_out->Nt = _Nt;
	lattice_out->totalV = _N*_N*_N*_Nt;
	lattice_out->N2 = _N*_N;
	lattice_out->spatialV = _N*_N*_N;
	lattice_out->dimLattice = totalV*4*4;
	lattice_out->colorV = 3*totalV;
	lattice_out->U = malloc(sizeof(double)*lattice_out->dimLattice);
  	return lattice_out;
}

LatticeColorVectorReal* newLatticeColorVectorReal(int _N, int _Nt) {    // emulates a constructor for lattice object 
    LatticeColorVectorReal* vector_out = malloc(sizeof(LatticeColorVectorReal));
    vector_out->N = _N;
    vector_out->Nt = _Nt;
	vector_out->dimVector = _N*_N*_N*_Nt*3;
    vector_out->vec = malloc(sizeof(double)*vector_out->dimVector);
  return vector_out;
}

LatticeColorVectorComplex* newLatticeColorVectorComplex(int _N, int _Nt) {    // emulates a constructor for lattice object 
    LatticeColorVectorComplex* vector_out = malloc(sizeof(LatticeColorVectorComplex));
    vector_out->N = _N;
    vector_out->Nt = _Nt;
	vector_out->dimVector = _N*_N*_N*_Nt*3;
    vector_out->vec = malloc(sizeof(double complex)*vector_out->dimVector);
  return vector_out;
}

LatticeGaugeSU2* newLatticeGaugeSU2(int _N, int _Nt) {    // emulates a constructor for lattice object 
	LatticeGaugeSU2 * output = malloc(sizeof(LatticeGaugeSU2));
	output->N = _N;
	output->N2 = _N*_N;
	output->Nt = _Nt;
	output->spatialV = _N*_N*_N;
	output->totalV = spatialV*_Nt;
	output->dimg = output->totalV*4;
	output->g = malloc(sizeof(double)*output->dimg);
	return output;
}

void freeLatticeLinkSU2(LatticeLinkSU2 * l) {    // emulates a constructor for lattice object 
    free(l->U);
	free(l);
}

void freeLatticeColorVectorReal(LatticeColorVectorReal * v) {    // emulates a constructor for lattice object 
    free(v->vec);
	free(v);
}

void freeLatticeColorVectorComplex(LatticeColorVectorComplex * v) {    // emulates a constructor for lattice object 
    free(v->vec);
	free(v);
}

void freeLatticeGaugeSU2(LatticeGaugeSU2 * l) {    // emulates a constructor for lattice object 
    free(l->g);
	free(l);
}


