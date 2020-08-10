#include "statistics.h"
#include "utilities.h"
#include "algebra.h"

#include "nrutil.h"

double naive_avg(double * dataSet, int dim){
    double sum = 0e0;
    int i;
    for(i=0 ; i<dim ; i++)
        sum += dataSet[i];
    return(sum/(double)(dim));
}

double naive_sampleVar(double * dataSet, int dim){
    double sum = 0e0;
    double avg = naive_avg(dataSet,dim);
    int i;
    for(i=0;i<dim;i++)
        sum += (dataSet[i]-avg)*(dataSet[i]-avg);
    return(sum/(double)dim);
}

int sortInt(int i , int j){
    //sort int in the interval [i,j)
    return((int)(ran0(global_seed)*(j-i) + i));
}

void plaquetteHistogram(double * histo, LatticeLinkSU2 * lattice, int nbins){
    //return histogram of plaquette values (-1,1) in a given lattice configuration

    int ibin;
    double start,final;
    double binwidth = 2.0/nbins;
    double sum = 0;
    //for(ibin=0;ibin<nbins;ibin++) histo[ibin] = 0;
    int temp;
    int k,j,mi,ni,t,x,y,z;
    int ani[4],ami[4];
    double * auxc = malloc(sizeof(double)*4);
    double * aux = malloc(sizeof(double)*4);
    int X=1;
    int Y=1;
    for(t=0 ; t<Nt ; t++){
        for(x=0 ; x<N ; x++){
            for(y=0 ; y<N ;y++){
                for(z=0 ; z<N ; z++){
                    for(ni=0 ; ni<3 ; ni++){
                        setUnitVector(&ani[0],ni);
                        for(mi=ni+1 ; mi<4 ; mi++){
                            setUnitVector(&ami[0],mi);
                            k=0;
                            j=0;
                            setidv(aux);

                            for(k=0 ; k<X ; k++){
								mmprodv(aux, aux, getLink(lattice, (t+k*ami[0])%Nt , (x+k*ami[1])%N , (y+k*ami[2])%N , (z+k*ami[3])%N ,mi) );
							}
							k--;

                            for(j=0 ; j<Y ; j++){
								mmprodv(aux , aux , getLink(lattice, (t+(k+1)*ami[0] + j*ani[0])%Nt , (x+(k+1)*ami[1]+ j*ani[1])%N ,
										(y+(k+1)*ami[2]+ j*ani[2])%N , (z+(k+1)*ami[3]+ j*ani[3])%N , ni) );
							}
							j--;

                            for( ; k > -1 ; k--){
								hermcv( auxc , getLink(lattice, (t+k*ami[0] + (j+1)*ani[0])%Nt , (x+k*ami[1]+ (j+1)*ani[1])%N ,
									(y+k*ami[2]+ (j+1)*ani[2])%N , (z+k*ami[3]+ (j+1)*ani[3])%N , mi ) );
								mmprodv(aux, aux , auxc);
							}
							k++;

							for( ; j > -1 ; j--){
								hermcv( auxc , getLink(lattice, (t+ j*ani[0])%Nt , (x+ j*ani[1])%N ,
									(y+ j*ani[2])%N , (z+ j*ani[3])%N , ni ) );
								mmprodv(aux , aux , auxc);
							}

                            //sum += tracev(aux)/2.;

                            for(ibin=0;ibin<nbins;ibin++){
                                start = -1 + ibin*binwidth;
                                final = -1 + (ibin+1)*binwidth;
                                if( (start <= tracev(aux)/2.) && (tracev(aux)/2. <= final) ){
                                    histo[ibin]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void bootstrapSimple(double * output, double * data, int n, int m, int K){
    double * resampled_avg = malloc(sizeof(double)*K);
    setzerovr(resampled_avg , K);
    int i,j;

    double bootstrap_avg = 0e0, bootstrap_variance = 0e0;
    for(i=0;i<K;i++){
        for(j=0;j<m;j++){
            resampled_avg[i] += data[sortInt(0,n)];
        }
        resampled_avg[i] /= (double)m;
        bootstrap_avg += resampled_avg[i];
    }
    bootstrap_avg /= (double)K;

    for(i=0;i<K;i++){
        bootstrap_variance += pow(resampled_avg[i]-bootstrap_avg,2);
    }
    bootstrap_variance /= (double)K;

    output[0] = bootstrap_avg;
    output[1] = sqrt(bootstrap_variance);       //converting to standard deviation

    //printf("\navgD = %lf +- %lf\n",output[0],output[1]);
    ////getchar();

    free(resampled_avg);
}

void bootstrapCreutzRatio(double * output, double * wij, double * wiMjM, double * wiMj, double * wijM, int n, int m, int K){
    //bootstrap average and std computation of a derived quantity chi(i,j) = -log(wij*wiMjM/wiMj*wijM)
    int i,j, index;
    double chi_bootstrap_avg = 0e0, chi_bootstrap_variance = 0e0;
    double * resampled_avg_wij = malloc(sizeof(double)*K);
    double * resampled_avg_wiMjM = malloc(sizeof(double)*K);
    double * resampled_avg_wiMj = malloc(sizeof(double)*K);
    double * resampled_avg_wijM = malloc(sizeof(double)*K);
    setzerovr(resampled_avg_wij , K);
    setzerovr(resampled_avg_wiMjM , K);
    setzerovr(resampled_avg_wiMj , K);
    setzerovr(resampled_avg_wijM , K);
    double * chi_resampled_avg = malloc(sizeof(double)*K);
    setzerovr(chi_resampled_avg , K);

    for(i=0;i<n;i++) printf("\n w = %lf",wij[i]);

    for(i=0;i<K;i++){
        for(j=0;j<m;j++){
            index = sortInt(0,n);
            resampled_avg_wij[i] += wij[index];
            resampled_avg_wiMjM[i] += wiMjM[index];
            resampled_avg_wiMj[i] += wiMj[index];
            resampled_avg_wijM[i] += wijM[index];
        }
        //printf("\n rw = %lf",resampled_avg_wij[i]);
        ////getchar();

        resampled_avg_wij[i] /= (double)m;
        resampled_avg_wiMjM[i] /= (double)m;
        resampled_avg_wiMj[i] /= (double)m;
        resampled_avg_wijM[i] /= (double)m;

        chi_resampled_avg[i] = -log( fabs(resampled_avg_wij[i]*resampled_avg_wiMjM[i]/(resampled_avg_wiMj[i]*resampled_avg_wijM[i]) ));
        chi_bootstrap_avg += chi_resampled_avg[i];
        if(isnan(chi_bootstrap_avg) || isnan(chi_resampled_avg[i])){
            printf("\n%lf %lf", chi_bootstrap_avg, chi_resampled_avg[i]);
            printf("\n%lf , %lf , %lf , %lf",resampled_avg_wij[i] , resampled_avg_wiMjM[i] , resampled_avg_wiMj[i] , resampled_avg_wijM[i]);
            ////getchar();
        }
    }
    chi_bootstrap_avg /= (double)K;

    for(i=0;i<K;i++){
        chi_bootstrap_variance += pow(chi_resampled_avg[i]-chi_bootstrap_avg,2);
        if(isnan(pow(chi_resampled_avg[i]-chi_bootstrap_avg,2)) ){
            printf("\n ih rapaz: %15lf %lf",chi_resampled_avg[i],chi_bootstrap_avg);
            printf("\n%15lf , %lf , %lf , %lf",resampled_avg_wij[i] , resampled_avg_wiMjM[i] , resampled_avg_wiMj[i] , resampled_avg_wijM[i]);
            ////getchar();
        }
    }
    chi_bootstrap_variance /= (double)K;

    output[0] = -log( fabs(naive_avg(wij,n)*naive_avg(wiMjM,n)/(naive_avg(wiMj,n)*naive_avg(wijM,n))) );
    output[1] = sqrt(chi_bootstrap_variance);       //converting to standard deviation

    if(isnan(output[1])){
        printf("\nNaN DUMP\n");
        for(i=0;i<n;i++){
            printf("\n%lf", chi_bootstrap_variance);
        }
        ////getchar();
    }

    printf("\navgD = %lf +- %lf\n",output[0],output[1]);

    free(chi_resampled_avg);
    free(resampled_avg_wij);
    free(resampled_avg_wiMjM);
    free(resampled_avg_wiMj);
    free(resampled_avg_wijM);
}

void bootstrapSimpleRatio(double * output, double * wij, double * wijM, int n, int m, int K){
    //bootstrap average and std computation of a derived quantity chi(i,j) = -log(wij*wiMjM/wiMj*wijM)
    int i,j, index;
    double chi_bootstrap_avg = 0e0, chi_bootstrap_variance = 0e0;
    double * resampled_avg_wij = malloc(sizeof(double)*K);
    double * resampled_avg_wijM = malloc(sizeof(double)*K);
    setzerovr(resampled_avg_wij , K);
    setzerovr(resampled_avg_wijM , K);
    double * chi_resampled_avg = malloc(sizeof(double)*K);
    setzerovr(chi_resampled_avg , K);

    for(i=0;i<K;i++){
        for(j=0;j<m;j++){
            index = sortInt(0,n);
            resampled_avg_wij[i] += wij[index];
            resampled_avg_wijM[i] += wijM[index];
        }
        resampled_avg_wij[i] /= (double)m;
        resampled_avg_wijM[i] /= (double)m;

        chi_resampled_avg[i] = -log( fabs(resampled_avg_wij[i]/resampled_avg_wijM[i]) );
        chi_bootstrap_avg += chi_resampled_avg[i];
        if(isnan(chi_bootstrap_avg)){
            printf("\n%lf", chi_bootstrap_avg);
            printf("\n%lf , %lf",resampled_avg_wij[i], resampled_avg_wijM[i]);
            ////getchar();
        }
    }
    chi_bootstrap_avg /= (double)K;

    for(i=0;i<K;i++){
        chi_bootstrap_variance += pow(chi_resampled_avg[i]-chi_bootstrap_avg,2);

    }
    chi_bootstrap_variance /= (double)K;

    output[0] = -log( fabs(naive_avg(wij,n)/naive_avg(wijM,n)) );
    output[1] = sqrt(chi_bootstrap_variance);       //converting to standard deviation

    if(isnan(output[1])){
        printf("\nNaN DUMP\n");
        for(i=0;i<n;i++){
            //printf("\n%lf", chi_bootstrap_variance);
        }
        ////getchar();
    }

    printf("\navgD = %lf +- %lf\n",output[0],output[1]);

    free(chi_resampled_avg);
    free(resampled_avg_wij);
    free(resampled_avg_wijM);
}

void bootstrapSimpleLog(double * output, double * wij, int n, int m, int K){
    //bootstrap average and std computation of a derived quantity gamma(i,j) = -log(W(i,j))
    int i,j, index;
    double gamma_bootstrap_avg = 0e0, gamma_bootstrap_variance = 0e0;
    double * resampled_avg_wij = malloc(sizeof(double)*K);
    setzerovr(resampled_avg_wij , K);
    double * gamma_resampled_avg = malloc(sizeof(double)*K);

    for(i=0;i<K;i++){
        for(j=0;j<m;j++){
            index = sortInt(0,n);
            resampled_avg_wij[i] += wij[index];
        }
        resampled_avg_wij[i] /= (double)m;

        gamma_resampled_avg[i] = -log(fabs(resampled_avg_wij[i]));
        gamma_bootstrap_avg += gamma_resampled_avg[i];
        if(isnan(gamma_bootstrap_avg)){
            printf("\nFOUND A BOOTSTRAP NaN = %lf", gamma_bootstrap_avg);
            printf("\n%lf ",resampled_avg_wij[i]);
            ////getchar();
        }
    }
    gamma_bootstrap_avg /= (double)K;

    for(i=0;i<K;i++){
        gamma_bootstrap_variance += pow(gamma_resampled_avg[i]-gamma_bootstrap_avg,2);
        if(isnan(gamma_bootstrap_variance)){
            //printf("DUMP: %lf | %lf ",gamma_resampled_avg[i],gamma_bootstrap_avg);
            ////getchar();
        }
    }
    gamma_bootstrap_variance /= (double)K;

    output[0] = -log(fabs(naive_avg(wij,n)));
    output[1] = sqrt(gamma_bootstrap_variance);       //converting to standard deviation



    free(gamma_resampled_avg);
    free(resampled_avg_wij);
}

void cornellf(double x, double f[], int np){
    f[1] = 1.0;
    f[2] = 1.0/x;
    f[3] = x;
}

//from Numerical Recipes 90: functions leading to linear fit()
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double gammln(double xx){
    //Returns the value ln[Γ( xx )] for xx > 0.
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void gcf(double *gammcf, double a, double x, double *gln){
    //Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction represen-
    //tation as gammcf . Also returns ln Γ(a) as gln .
    double gammln(double xx);
    void nrerror(char error_text[]);
    int i;
    double an,b,c,d,del,h;
    *gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) {
        an = -i*(i-a);
        b += 2.0;
        d=an*d+b;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=b+an/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
    *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(double *gamser, double a, double x, double *gln){
    //Returns the incomplete gamma function P (a, x) evaluated by its series representation as gamser .
    //Also returns ln Γ(a) as gln .
    double gammln(double xx);
    void nrerror(char error_text[]);
    int n;
    double sum,del,ap;
    *gln=gammln(a);
    if (x <= 0.0){
        if (x < 0.0) nrerror("x less than 0 in routine gser");
        *gamser=0.0;
        return;
    }
    else{
        ap=a;
        del=sum=1.0/a;
        for (n=1;n<=ITMAX;n++){
            ++ap;
            del *= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS){
                *gamser=sum*exp(-x+a*log(x)-(*gln));
                return;
            }
        }
        nrerror("a too large, ITMAX too small in routine gser");
        return;
    }
}

double gammq(double a, double x){
    //Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
    void gcf(double *gammcf, double a, double x, double *gln);
    void gser(double *gamser, double a, double x, double *gln);
    void nrerror(char error_text[]);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
    if(x < (a+1.0)){
        gser(&gamser,a,x,&gln);
        return 1.0-gamser;
    }
    else{
        gcf(&gammcf,a,x,&gln);
        return gammcf;
    }
}

void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
    double *b, double *siga, double *sigb, double *chi2, double *q){
    /*
    Given a set of data points x[1..ndata] , y[1..ndata] with individual standard deviations
    sig[1..ndata] , fit them to a straight line y = a + bx by minimizing χ 2 . Returned are
    a,b and their respective probable uncertainties siga and sigb , the chi-square chi2 , and the
    goodness-of-fit probability q (that the fit would have χ 2 this large or larger). If mwt=0 on
    input, then the standard deviations are assumed to be unavailable: q is returned as 1.0 and
    the normalization of chi2 is to unit standard deviation on all points.
    */
    double gammq(double a, double x);
    int i;

    double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
    *b=0.0;
    if(mwt){
        ss=0.0;
        for(i=0;i<ndata;i++){
            wt=1.0/SQR(sig[i]);
            ss += wt;
            sx += x[i]*wt;
            sy += y[i]*wt;
        }
    }
    else{
        for(i=0;i<ndata;i++){
            sx += x[i];
            sy += y[i];
        }
        ss=ndata;
    }
    sxoss=sx/ss;
    if(mwt){
        for(i=0;i<ndata;i++){
            t=(x[i]-sxoss)/sig[i];
            st2 += t*t;
            *b += t*y[i]/sig[i];
        }
    }
    else{
        for(i=0;i<ndata;i++){
            t=x[i]-sxoss;
            st2 += t*t;
            *b += t*y[i];
        }
    }
    *b /= st2;
    *a=(sy-sx*(*b))/ss;
    *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    *sigb=sqrt(1.0/st2);
    *chi2=0.0;
    *q=1.0;
    if(mwt == 0){
        for(i=0;i<ndata;i++)
            *chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
        sigdat=sqrt((*chi2)/(ndata-2));
        *siga *= sigdat;
        *sigb *= sigdat;
    }
    else{
        for(i=0;i<ndata;i++)
            *chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
        if(ndata>2)
            *q=gammq(0.5*(ndata-2),0.5*(*chi2));

    }
}

void svdvar(double **v, int ma, double w[], double **cvm){
    //To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma parameters obtained
    //by svdfit , call this routine with matrices v[1..ma][1..ma] , w[1..ma] as returned from
    //svdfit.
    int k,j,i;
    double sum,*wti;
    wti=vector(1,ma);
    for (i=1;i<=ma;i++) {
        wti[i]=0.0;
        if (w[i]) wti[i]=1.0/(w[i]*w[i]);
    }
    for (i=1;i<=ma;i++) {
        for (j=1;j<=i;j++) {
            for (sum=0.0,k=1;k<=ma;k++) sum += v[i][k]*v[j][k]*wti[k];
            cvm[j][i]=cvm[i][j]=sum;
        }
    }
    free_vector(wti,1,ma);
}

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]){
    //Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n] , w[1..n] ,
    //v[1..n][1..n] as returned by svdcmp . m and n are the dimensions of a , and will be equal for
    //square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
    //No input quantities are destroyed, so the routine may be called sequentially with different b ’s.
    int jj,j,i;
    double s,*tmp;
    tmp=vector(1,n);
    for (j=1;j<=n;j++) {
        s=0.0;
        if (w[j]) {
            for (i=1;i<=m;i++) s += u[i][j]*b[i];
            s /= w[j];
        }
        tmp[j]=s;
    }
    for (j=1;j<=n;j++) {
        s=0.0;
        for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
        x[j]=s;
    }
    free_vector(tmp,1,n);
}

void svdcmp(double **a, int m, int n, double w[], double **v){
    //Given a matrix a[1..m][1..n] , this routine computes its singular value decomposition, A =
    //U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is out-
    //put as a vector w[1..n] . The matrix V (not the transpose V T ) is output as v[1..n][1..n] .
    double pythag(double a, double b);
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
    rv1=vector(1,n);
    g=scale=anorm=0.0;
    for (i=1;i<=n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m) {
            for (k=i;k<=m;k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k=i;k<=m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<=m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i <= m && i != n) {
            for (k=l;k<=n;k++) scale += fabs(a[i][k]);
            if (scale) {
                for (k=l;k<=n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<=m;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l;k<=n;k++) a[i][k] *= scale;
            }
        }
        anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) {
        if (i < n) {
            if (g) {
                for (j=l;j<=n;j++)
                v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=IMIN(m,n);i>=1;i--) {
        l=i+1;
        g=w[i];
        for (j=l;j<=n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            for (j=l;j<=n;j++) {
                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<=m;j++) a[j][i] *= g;
        } else for (j=i;j<=m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n;k>=1;k--) {
        for (its=1;its<=30;its++) {
            flag=1;
            for (l=k;l>=1;l--) {
                nm=l-1;
                if ((double)(fabs(rv1[l])+anorm) == anorm) {
                    flag=0;
                    break;
                }
                if ((double)(fabs(w[nm])+anorm) == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((double)(fabs(f)+anorm) == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=1;j<=m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=1;jj<=n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=1;jj<=m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    free_vector(rv1,1,n);
}

double pythag(double a, double b){
    //Computes (a 2 + b 2 ) 1/2 without destructive underflow or overflow.
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

#define TOL 1.0e-5

void svdfit(double x[], double y[], double sig[], int ndata, double a[],
     int ma, double **u, double **v, double w[], double *chisq,
        void (*funcs)(double, double [], int)){
    //Given a set of data points x[1..ndata] , y[1..ndata] with individual standard deviations
    //sig[1..ndata] , use χ minimization to determine the coefficients a[1..ma] of the fitting
    //function y = ai i × afunc i (x). Here we solve the fitting equations using singular
    //value decomposition of the ndata by ma matrix, as in §2.6. Arrays u[1..ndata][1..ma] ,
    //v[1..ma][1..ma] , and w[1..ma] provide workspace on input; on output they define the
    //singular value decomposition, and can be used to obtain the covariance matrix. The pro-
    //gram returns values for the ma fit parameters a , and χ 2 , chisq . The user supplies a routine
    //funcs(x,afunc,ma) that returns the ma basis functions evaluated at x = x in the array afunc[1..ma] .
    int j,i;
    void svbksb(double **u, double w[], double **v, int m, int n, double b[],
    double x[]);
    void svdcmp(double **a, int m, int n, double w[], double **v);
    double wmax,tmp,thresh,sum,*b,*afunc;
    b=vector(1,ndata);
    afunc=vector(1,ma);
    for (i=1;i<=ndata;i++) {
        (*funcs)(x[i],afunc,ma);
        tmp=1.0/sig[i];
        for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
        b[i]=y[i]*tmp;
    }
    svdcmp(u,ndata,ma,w,v);
    wmax=0.0;
    for (j=1;j<=ma;j++)
    if (w[j] > wmax) wmax=w[j];
    thresh=TOL*wmax;
    for (j=1;j<=ma;j++)
    if (w[j] < thresh) w[j]=0.0;
    svbksb(u,w,v,ndata,ma,b,a);
    *chisq=0.0;
    for (i=1;i<=ndata;i++) {
        (*funcs)(x[i],afunc,ma);
        for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
        *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
    }
    free_vector(afunc,1,ma);
    free_vector(b,1,ndata);
}

//data blocking
void dataBlockingOnlyDivisors(double * data, int dim){
    FILE * b;
    int i, binsize, nbins, ibin;
    int divisorsCount=0;

    //calculating the number of divisors=binvar array size
    for(binsize=1;binsize<(int)(dim/2.0);binsize++){
        if(dim%binsize==0)
            divisorsCount++;
    }
    double * binsizeList = malloc(sizeof(double)*divisorsCount);
    double * binvar = malloc(sizeof(double)*divisorsCount);

    //executing blocking
    double * bindata;
    ibin=0;
    for(binsize=1;binsize<(int)(dim/2.0);binsize++){
        nbins = (int)(dim/(double)binsize);
        if(dim%binsize==0){
            bindata = malloc(sizeof(double)*nbins);
            //printf("\n binsize=%d ; nbins=%d",binsize,nbins);
            for(i=0;i<nbins;i++){                                               //creating binned data
                bindata[i] = naive_avg(data+i*binsize,binsize);
                //printf("\n - i=%d - avg=%lf",i,naive_avg((data+i),binsize));
            }
            binvar[ibin] = naive_sampleVar(bindata,nbins);
            binsizeList[ibin] = binsize;
            printf("\n%lf %lf ",binsizeList[ibin],binvar[ibin]);
            //fprintf(b,"\n%e %e",binsizeList[ibin],binvar[ibin]);
            fflush(stdout);
            ibin++;
            free(bindata);
        }
    }
    b = fopen("./data/binvar.dat", "w");
    for(ibin=1;ibin<=divisorsCount;ibin++){
        fprintf(b,"\n%e %e",binsizeList[divisorsCount-ibin],binvar[divisorsCount-ibin]);
    }
    fclose(b);

    /*

    //quality of fitting to 1/x
    void blockingFitFunction(double x, double f[], int np){
        f[1] = 1.0;
        f[2] = 1.0/x;
    }
    double * sigmaDummy = malloc(sizeof(double)*divisorsCount);
    for(i=0;i<divisorsCount;i++) sigmaDummy[i]=0.1;
    double * fitParameters = vector(1,2);
    double **covarM = dmatrix(1,2,1,2);
    double **u;
    double **v = dmatrix(1,2,1,2);
    double  *w = vector(1,2);
    double chi2;
    int tempDataSize;
    i=0;
    for(i=0;i<divisorsCount;i++){
        tempDataSize = divisorsCount-i;

        u = dmatrix(1,tempDataSize,1,2);
        //svdfit(binsizeList, binvar, sigmaDummy, tempDataSize,
            //fitParameters, 2, u, v, w, &chi2, blockingFitFunction);
        //svdvar(v, 2, w, covarM);                                                //calculate covariance matrix
        //free_dmatrix(u,1,tempDataSize,1,2);

        //printf("\n popping %3d elements: fit(x)=%lf + %lf/x",i,fitParameters[1],fitParameters[2]);

    }
    */
    free(binvar);

}
