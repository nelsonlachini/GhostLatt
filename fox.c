#include "fox.h"
#include "global.h"

double xymapFox(double y , double x0){
	y -= x0;
	return( N*N*N*(1e0/((N-y)*(N-y)) - 1e0/(y*y)) );
}

double xymapRadius(double y , double x0){
	return(y);
	//y += x0;
	//return( xymapFox(y + N/2e0 , x0) );
}

int isoutside(double t,double x,double y,double z, double x0, int M){
	t-=x0;
	x-=x0;
	y-=x0;
	z-=x0;
	double sup = N/2e0 + M, inf = N/2e0 - M;
	if( (t<=inf||t>=sup)||(x<=inf||x>=sup)||(y<=inf||y>=sup)||(z<=inf||z>=sup) ){
		return(1);
	}
	else{
		return(0);
	}
}

void amuregVFox(double * amu, double t, double x, double y, double z, int mu, double x0,double R){
	double xv[4] = {xymapFox(t,x0),xymapFox(x,x0),xymapFox(y,x0),xymapFox(z,x0)};
	double x2 = xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2] + xv[3]*xv[3];
	double rho = xymapRadius(R,x0);

	amu[0] = 0e0;

	if(mu==0){
		amu[1] = -xv[1];
		amu[2] = -xv[2];
		amu[3] = -xv[3];
	}
	else if(mu==1){
		amu[1] = xv[0];
		amu[2] = xv[3];
		amu[3] = -xv[2];
	}
	else if(mu==2){
		amu[1] = -xv[3];
		amu[2] = xv[0];
		amu[3] = xv[1];
	}
	else if(mu==3){
		amu[1] = xv[2];
		amu[2] = -xv[1];
		amu[3] = xv[0];
	}

	cmprodv(amu,1e0/(x2+rho*rho),amu);
}

void amusingVFox(double * amu, double t , double x, double y, double z, int mi, double x0, double R){
	int i,k;
	double xv[4] = {xymapFox(t,x0),xymapFox(x,x0),xymapFox(y,x0),xymapFox(z,x0)};
	double * sum = malloc(sizeof(double)*4);
	double C = xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2] + xv[3]*xv[3];

	double rho = xymapRadius(R,x0);
	C = rho*rho/(C*(C+rho*rho));

	//if(isnan(xymapRadius(R,x0))){
	//	printf("\n %lf",xymapRadius(R,x0));
//		getchar();
	//}





	if(t-x0==0 || x-x0==0 || y-x0==0 || z-x0==0){
		setzerov(amu);
	}
	else{
		if(mi==0){
			sum[0] = 0e0;
			for(i=1;i<4;i++){
				sum[i] = xv[i];
			}
		}
		else{
			setzerov(sum);
			sum[0] = xv[0];
			for(i=1;i<4;i++){
				for(k=1;k<4;k++){
					sum[k] += eps4(0,i,mi,k)*xv[i];
				}
			}
		}
		cmprodv(amu,C,sum);
	}
	free(sum);
}

void umuregVFox(double * umu, int t, int x, int y, int z, int mi, double x0, double R){
	double anorm;
	int i;
	double * aux = malloc(sizeof(double)*4);
	double * amu = malloc(sizeof(double)*4);
	double emi[4] = {0,0,0,0};
	int i_number = 1;
	double i_step = 1e0/(i_number+1);

	setzerov(amu);
	for(i=0;i<i_number;i++){
		emi[mi] += i_step;
		amuregVFox(aux,getStep(t,emi[0]),getStep(x,emi[1]),getStep(y,emi[2]),getStep(z,emi[3]),mi,x0,R);
		sumv(amu,amu,aux);
	}
	cmprodv(amu,1e0/i_number,amu);

	anorm = sqrt(amu[1]*amu[1] + amu[2]*amu[2] + amu[3]*amu[3]);
	if(anorm!=0){
		umu[0] = cos(anorm);
		for(i=1;i<4;i++){
			umu[i] = amu[i]*sin(anorm)/anorm;
		}
	}
	else{
		setidv(umu);
	}

	free(aux);
	free(amu);
}

void umusingVFox(double * umu, int t, int x, int y, int z, int mi, double x0, double R){
	double anorm;
	int i,ami[4]={0,0,0,0};
	double * aux = malloc(sizeof(double)*4);
	double * amu = malloc(sizeof(double)*4);

	setzerov(amu);

	//amusingVFox(aux,t,x,y,z,mi);
	//sumv(amu,aux,amu);

	ami[mi] = 1;
	amusingVFox(aux,getStep(t,ami[0]/2e0),getStep(x,ami[1]/2e0),getStep(y,ami[2]/2e0),getStep(z,ami[3]/2e0),mi,x0,R);

	sumv(amu,aux,amu);

	//amusingVFox(aux,getStep(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]),mi);
	//sumv(amu,aux,amu);

	//cmprodv(amu,1e0/3e0,amu);
	//cmprodv(amu,5e-1,amu);



	anorm = sqrt(amu[1]*amu[1] + amu[2]*amu[2] + amu[3]*amu[3]);



	if(anorm!=0){
		umu[0] = cos(anorm);
		for(i=1;i<4;i++){
			umu[i] = amu[i]*sin(anorm)/anorm;
		}
	}
	else{
		setidv(umu);
	}



	free(aux);
	free(amu);
}

void creategFoxV(double * g, double t, double x, double y, double z, double x0){
	double * aux = malloc(sizeof(double)*4);
	double xv[4] = {xymapFox(t,x0),xymapFox(x,x0),xymapFox(y,x0),xymapFox(z,x0)};
	double xnorm = sqrt(xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2] + xv[3]*xv[3]);
	int i;

	for(i=0;i<4;i++){
		g[i] = xv[i];
	}
	cmprodv(g,1e0/xnorm,g);
	free(aux);
}

void foxInstanton(double * U , double R, double x0 , int M){
	//using umusing
	int t,x,y,z,mu,i;
	int ami[4];
	double * auxg = malloc(sizeof(double)*4);
	double * auxu = malloc(sizeof(double)*4);
	//int c=0;
	for(t=0 ; t<N ; t++){
		for(x=0 ; x<N ; x++){
			for(y=0 ; y<N ; y++){
				for(z=0 ; z<N ; z++){
					for(mu=0 ; mu<4 ; mu++){
						setquadv(&ami[0],mu);

						if( isoutside(t,x,y,z,x0,M) && isoutside(getStep(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]),x0,M) ){
							//outside
							umusingVFox(auxu,t,x,y,z,mu,x0,R);
							copyv( getU(U,t,x,y,z,mu) , auxu);
						}
						else{
							//inside
							umuregVFox(auxu,t,x,y,z,mu,x0,R);

							if( isoutside(t,x,y,z,x0,M) ){
								//origin of the link is in the interface
								creategFoxV(auxg,t,x,y,z,x0);
								hermcv(auxg,auxg);
								mmprodv(auxu,auxg,auxu);
							}
							else if( isoutside(getStep(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]),x0,M) ){
								//endpoint is in the interface
								creategFoxV(auxg,getStep(t,ami[0]),getStep(x,ami[1]),getStep(y,ami[2]),getStep(z,ami[3]),x0);
								mmprodv(auxu,auxu,auxg);
							}
							copyv( getU(U,t,x,y,z,mu) , auxu);
							//printf("\n Hey, I'm here!!");
							//printpos(t,x,y,z,mu);
							//printc(isoutside(0,0,1,0,x0,M));
							//getchar();
						}

						if(isnan(detv(getU(U,t,x,y,z,mu)))){
							printpos(t,x,y,z,mu);

							getchar();
						}
					}
				}
			}
		}
	}
	free(auxu);
	free(auxg);
	//printf("\n %d \n",c);
}
