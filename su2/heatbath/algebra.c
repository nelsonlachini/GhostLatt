
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

void cmprodv( double * pout , double complex k, double * C ){
	double A[4];
	double a = creal(k);
	double b = cimag(k);
	int i;
	
	A[0] = a*(*getelv(C,0)) - b*(*getelv(C,3));
	A[1] = a*(*getelv(C,1)) - b*(*getelv(C,2));
	A[2] = a*(*getelv(C,2)) + b*(*getelv(C,1));	
	A[3] = a*(*getelv(C,3)) + b*(*getelv(C,0));
	
	for(i=0 ; i<4 ; i++){
		*(getelv(pout,i)) = A[i];
	}
}

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

void setzerov(double * vector){
	int i;
	for (i=0 ; i<4 ; i++){
		*getelv(vector , i) = 0.0;
	}
}

int signal(double number){
	//returns +1 if element is (zero) positive and -1 else
	if(number < 0)	return -1;	
	return 1;
}
