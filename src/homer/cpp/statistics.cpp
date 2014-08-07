/*
 * This is a file devoted to statistics
 */
#include "statistics.h"

double logPoisson(int x, double lambda) {
	double rv = ((double)x)*log(lambda) - lambda - factln(x);
	return rv;
}
double logCumulativePoisson(int x, double lambda) {
	double P = 0.0;
	for (int i=0;i<=x;i++) {
		double p = logPoisson(i, lambda);
		if (i==0) {
			P = p;
		} else {
			P += log(1+exp(p-P));
		}
	}
	return P;
}
double ilogCumulativePoisson(int x, double lambda) {
	double P = 0.0;
	double diff = 1;
	int limit = (int)lambda*10+x;
	for (int i=x;i<=limit;i++) {
		double p = logPoisson(i, lambda);
		if (i==x) {
			P = p;
		} else {
			diff = log(1+exp(p-P));
			P += diff;
		}
		if (diff <= 0) break;
	}
	return P;
}


//threshold helps us determine when to stop bothering
double cumulativePoisson(int x, double lambda) {
	double sum = 0.0;
	for (int i=0;i<=x;i++) {
		sum += exp(logPoisson(i, lambda));
	}
	return sum;
}


StatMemory::StatMemory(int size, int NN, int nn1) {
	maxBytes = size;
	N = NN;
	n1 = nn1;
	denominator = bicoln(N,n1);
	length = (int)floor(sqrt((double)size)*2)-1;
	curBytes = 0;
	result = new float*[length];
	for (int i=0;i<length;i++) {
		result[i] = NULL;
	}
}

StatMemory::~StatMemory() {
	for (int i=0;i<length;i++) {
		if (result[i] != NULL) {
			delete [](result[i]);
		}
	}
	delete []result;
}

float StatMemory::getStat(int n2, int n) {
	// needs to grow if it's not big enough
	if (n2+1 > length) {
		float ** newresult = new float*[n2+1];
		for (int i=0;i<length;i++) {
			newresult[i] = result[i];
		}
		for (int i=length;i<n2+1;i++) {
			newresult[i] = NULL;
		}
		delete []result;
		result = newresult;
	}

	if (result[n2] != NULL) {
		float r = result[n2][n];
		if (r < EMPTY_STAT-1.0) return r;
	}

	if (result[n2] == NULL && curBytes < maxBytes) {
		int m = n2+1;
		if (m > n1+1) m = n1+1;
		result[n2] = new float[m];
		for (int i=0;i<m;i++) {
			result[n2][i] = EMPTY_STAT;
		}
		curBytes += m;
	}
	
    float P=0;
    int minIndex = (n1>n2)?n2:n1;

	for (int i=minIndex;i>=n;i--) {
	//for (int i=n;i<=minIndex;i++) { 

		if (result[n2] != NULL && result[n2][i] < EMPTY_STAT-1.0) {
			P = result[n2][i];
		} else {
        	float c1 = bicoln(n2,i);
        	float c2 = bicoln(N-n2,n1-i);
        	float p = ((c2+c1)-denominator);
			if (i==minIndex) {
				P = p;
			} else {
				P += log(1+exp(p-P));
			}
		}
		if (result[n2] != NULL && result[n2][i] > EMPTY_STAT-1.0) {
			result[n2][i] = P;
		}
    }
	//exit(0);
	return P;
}



float rankSumStat(int ranksum, int n1, int N) {
	if (n1 == 0) {
		return FLT_MAX;
	}

	float n2 = (float)(N-n1);
	float k1 = (float)ranksum - ((float)(n1*(n1+1.0)))/2.0;

	float z = ((float)((n1*n2)/2.0) - k1) 
					/ sqrt((n1*n2*(n1+n2+1.0))/12.0);
	return z;

/*
	float best = (float)(N*(N+1.0)/2.0);
	float ranksum2 = best - ranksum;
	float n2 = (float)(N - n1);
	//float stat = ((float)(ranksum*ranksum))/(float)n1;
	//if (((float)ranksum)/((float)n1) < best/((float)N)) {
	//	return stat;
	//} else {
	//	return -1*stat;
	//}

	float a1 = ((float)ranksum-(float)n1*((float)N+1.0)/2.0);
	float stat = a1*a1/(float)n1;
	float a2 = ((float)ranksum2-(float)n2*((float)N+1.0)/2.0);
	stat += a2*a2/n2;
	stat *= 12/(float)((N*(N+1)));
	float pvalue = chi2pvalue(stat, 1);
	if (((float)ranksum)/((float)n1) < best/((float)N)) {
		pvalue = 1.0 - pvalue;
	}
	return pvalue;
*/	
	
}

float ttest(float* data1, int n1, float* data2, int n2, float &t) {
    float var1,var2,svar,df,ave1,ave2,prob;

    avevar(data1,n1,ave1,var1);
    avevar(data2,n2,ave2,var2);

    df=n1+n2-2;
    svar=((n1-1)*var1+(n2-1)*var2)/df;
    t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
    prob=betai(0.5*df,0.5,df/(df+t*t));
    return prob;
}

void avevar(float* data, int n, float &ave, float &var) {
    int j;
    float s,ep;
    for (ave=0.0,j=0;j<n;j++) ave+= data[j];
    ave/=n;
    var = ep = 0.0;
    for (j=0;j<n;j++) {
        s = data[j]-ave;
        ep += s;
        var += s*s;
    }
    var = (var-ep*ep/n)/(n-1);
}

float logbetai(float a, float b, float x) {
    if (x <= 0.0 || x >= 1.0) {
        fprintf(stderr, "X is not right in betai (%f)\n",x);
        //exit(0);
    }
    //if (x==0.0 || x == 1.0) bt = -1*FLT_MAX;
    float bt = gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x);
    if (x < (a+1.0)/(a+b+2.0))
        return bt + log(betacf(a,b,x))-log(a);
    else
        return log(1.0-exp(bt)*betacf(b,a,1.0-x)/b);
}

double logbetaiD(double a, double b, double x) {
    if (x <= 0.0 || x >= 1.0) {
        fprintf(stderr, "X is not right in betai (%lf)\n",x);
        //exit(0);
    }
    //if (x==0.0 || x == 1.0) bt = -1*FLT_MAX;
    double bt = gammlnD(a+b)-gammlnD(a)-gammlnD(b)+a*log(x)+b*log(1.0-x);
    if (x < ((a+1.0)/(a+b+2.0))) {
        return bt + log(betacfD(a,b,x))-log(a);
	} else {
        return log(1.0-exp(bt)*betacfD(b,a,1.0-x)/b);
	}
}


float betai(float a, float b, float x) {
    float bt;
    if (x < 0.0 || x > 1.0) {
        fprintf(stderr, "X is not right in betai\n");
        //exit(0);
    }
    if (x==0.0 || x == 1.0) bt =0.0;
    else bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0))
        return bt*betacf(a,b,x)/a;
    else
        return 1.0-bt*betacf(b,a,1.0-x)/b;
}

#define MAXIT 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30

float betacf(float a, float b, float x) {
    int m,m2;
    float aa, c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h*=del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT) {
        fprintf(stderr, "Error calculating betacf\n");
        //exit(0);
    }
    return h;
}

//#define EPS 3.0e-7
#define FPMIND 1.0e-30
double betacfD(double a, double b, double x) {
    int m,m2;
    double aa, c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIND) d=FPMIND;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIND) d=FPMIND;
        c=1.0+aa/c;
        if (fabs(c) < FPMIND) c=FPMIND;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIND) d=FPMIND;
        c=1.0+aa/c;
        if (fabs(c) < FPMIND) c=FPMIND;
        d=1.0/d;
        del=d*c;
        h*=del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT) {
        fprintf(stderr, "Error calculating betacf\n");
        //exit(0);
    }
    return h;
}


float gammln(float xx) {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5396239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

//float bicoln(int n, int k);
//float factln(int n);
//float gammln(float xx);
//float hypergeo(int N, int n1, int n2, int n);

float hypergeo(int N, int n1, int n2, int n) {
    float P=0,c1=0,c2=0,c3=0;
    int i,minIndex;
    minIndex = (n1>n2)?n2:n1;
    c3 = bicoln(N,n2);
    for (i=n;i<=minIndex;i++) { 
        c1 = bicoln(n1,i);
        c2 = bicoln(N-n1,n2-i);
        P += exp((c2+c1)-c3);
    }
	return P;
}

float loghypergeo(int N, int n1, int n2, int n) {
//fprintf(stderr, "%d,%d,%d,%d\n", N,n1,n2,n);
    float P=0,c1=0,c2=0,c3=0;
    int i,minIndex;
    minIndex = (n1>n2)?n2:n1;
    c3 = bicoln(N,n2);
    for (i=n;i<=minIndex;i++) { 
        c1 = bicoln(n1,i);
        c2 = bicoln(N-n1,n2-i);
        float p = ((c2+c1)-c3);
		//printf("%e\n",p);
		if (i==n) {
			P = p;
		} else {
			P += log(1+exp(p-P));
		}
    }
	//exit(0);
	return P;
}


float iloghypergeo(int N, int n1, int n2, int n) {
//fprintf(stderr, "%d,%d,%d,%d\n", N,n1,n2,n);
    float P=0,c1=0,c2=0,c3=0;
	int lowerLimit = n1+n2-N;
	if (lowerLimit < 0) lowerLimit = 0;
    c3 = bicoln(N,n2);
    for (int i=n;i>=lowerLimit;i--) { 
        c1 = bicoln(n1,i);
        c2 = bicoln(N-n1,n2-i);
        float p = ((c2+c1)-c3);
		//printf("%e\n",p);
		if (i==n) {
			P = p;
		} else {
			P += log(1+exp(p-P));
		}
    }
	//exit(0);
	return P;
}


// double routines

double loghypergeoD(unsigned int N, unsigned int n1, unsigned int n2, unsigned int n) {
//fprintf(stderr, "%d,%d,%d,%d\n", N,n1,n2,n);
    double P=0,c1=0,c2=0,c3=0;
    int i,minIndex;
    minIndex = (n1>n2)?n2:n1;
    c3 = bicolnD(N,n2);
    for (i=n;i<=minIndex;i++) { 
        c1 = bicolnD(n1,i);
        c2 = bicolnD(N-n1,n2-i);
        double p = ((c2+c1)-c3);
		//printf("%e\n",p);
		if (i==(int)n) {
			P = p;
		} else {
			P += log(1+exp(p-P));
		}
    }
	return P;
}
double iloghypergeoD(unsigned int N, unsigned int n1, unsigned int n2, unsigned int n) {
//fprintf(stderr, "%d,%d,%d,%d\n", N,n1,n2,n);
    double P=0,c1=0,c2=0,c3=0;
    int i=0;
	int lowerLimit = n1+n2-N;
	if (lowerLimit < 0) lowerLimit = 0;
    c3 = bicolnD(N,n2);
    for (i=n;i>=lowerLimit;i--) { 
        c1 = bicolnD(n1,i);
        c2 = bicolnD(N-n1,n2-i);
        double p = ((c2+c1)-c3);
		//printf("%e\n",p);
		if (i==(int)n) {
			P = p;
		} else {
			P += log(1+exp(p-P));
		}
    }
	return P;
}


/*
float loghypergeoasdf(int N, int n1, int n2, int n) {
//fprintf(stderr, "%d,%d,%d,%d\n", N,n1,n2,n);
    float P=0,c1=0,c2=0,c3=0;
    int i,minIndex;
    minIndex = (n1>n2)?n2:n1;
    c3 = bicoln(N,n2);
	static float store[100];

	int setSize = 4;
	int maxAcc = setSize*4;
	if (minIndex - n < maxAcc) {
    	for (i=n;i<=minIndex;i++) { 
    	    c1 = bicoln(n1,i);
    	    c2 = bicoln(N-n1,n2-i);
    	    float p = ((c2+c1)-c3);
			//printf("%e\n",p);
			if (i==n) {
				P = p;
			} else {
				P += log(1+exp(p-P));
			}
    	}
	} else {
		int incFlag = 0;
		int decFlag = 0;
		int curIndex = 0;
		int retry = 0;
    	for (i=n;i<=minIndex;i++) { 
			retry++;
			int doneFlag = 0;
    	    c1 = bicoln(n1,i);
    	    c2 = bicoln(N-n1,n2-i);
    	    float p = ((c2+c1)-c3);
			store[curIndex] = p;
			if (i!=n) {
				float lastvalue = 0;
				if (curIndex == 0) {
					lastvalue = store[99];
				} else {
					lastvalue = store[curIndex-1];
				}
				if (p < lastvalue) {
					decFlag = 1;
				}
				//compute ratio
				if (i-n >= 2*setSize && retry >= setSize) {
					retry = 0;
					float baseline = p;
					float num=0.0,den=0.0;
					for (int i=0;i<2*setSize;i++) {
						int index = curIndex-i;
						if (index < 0) {
							index += 100;
						}
						if (i > setSize) {
							num += exp(store[index]-p);
						} else {
							den += exp(store[index]-p);
						}
					}
					float ratio = log(num/den);
					//fprintf(stdout, "%f\n", ratio);
					if (fabs(ratio) > 4) {
						if (decFlag) {
							doneFlag = 1;
						} else {
							P = 0.0;
							doneFlag = 1;
						}
					}
				}
			}
			curIndex++;
			if (curIndex >= 100) curIndex = 0;
			if (i==n) {
				P = p;
			} else {
				P += log(1+exp(p-P));
			}
			if (doneFlag == 1) {
				break;
			}
    	}
	}
	return P;
}
*/

float logbinomial(int n, int k, float r, int N) {
	float Llimit = ((float)1.0) /((float)N);
	float Hlimit = ((float)(N-1.0))/((float)N);
	if (r < Llimit) {
		r = Llimit;
	}
	if (r > Hlimit) {	
		r = Hlimit;
	}
	if (k == 0) {
		return 0.0;
	}
	return logbetai((float)k,(float)(n-k+1),(float)r);
}

double logbinomialD(unsigned int n, unsigned int k, double r, unsigned int N) {
	double Llimit = ((double)1.0) /((double)N);
	double Hlimit = ((double)(N-1.0))/((double)N);
	if (r < Llimit) {
		r = Llimit;
	}
	if (r > Hlimit) {	
		r = Hlimit;
	}
	if (k == 0) {
		return 0.0;
	}
	return logbetaiD((double)k,(double)(n-k+1),(double)r);
}

double ilogbinomialD(unsigned int n, unsigned int k, double r, unsigned int N) {
//fprintf(stderr, "input = %d, %d, %le, %d\n", n,k,r,N);
	double Llimit = ((double)1.0) /((double)N);
	double Hlimit = ((double)(N-1.0))/((double)N);
	if (r < Llimit) {
		r = Llimit;
	}
	if (r > Hlimit) {	
		r = Hlimit;
	}
	r = 1.0-r;
	k = n-k;

	if (k == 0) {
		return 0.0;
	}
	double rv = logbetaiD((double)k,(double)(n-k+1),(double)r);
//fprintf(stderr, "input = %d, %d, %le, %d == %lf\n", n,k,r,N,rv);
	return rv;
}


float bicoln(int nn,int k) {
    //float v = floor(0.5+exp(factln(nn)-factln(k)-factln(nn-k)));
	//if (nn < 0) fprintf(stderr, "nn is less than 0 %d\n", nn);
	//if (k < 0) fprintf(stderr, "k is less than 0 %d\n", k);
	//if (nn-k < 0) fprintf(stderr, "nn-k is less than 0 n=%d, k=%d, %d\n",nn,k, nn-k);
    float v = factln(nn)-factln(k)-factln(nn-k);
	return v;
}
double bicolnD(unsigned int nn,unsigned int k) {
    //float v = floor(0.5+exp(factln(nn)-factln(k)-factln(nn-k)));
	//if (nn < 0) fprintf(stderr, "nn is less than 0 %d\n", nn);
	//if (k < 0) fprintf(stderr, "k is less than 0 %d\n", k);
	//if (nn-k < 0) fprintf(stderr, "nn-k is less than 0 n=%d, k=%d, %d\n",nn,k, nn-k);
    return factlnD(nn)-factlnD(k)-factlnD(nn-k);
}

float factln(int nn) {
    static float a[HYPERGEO_CACHE_SIZE+1];
    if (nn<0) { 
		fprintf(stderr,"negative factorial %d\n", nn);
		//fprintf(stdout,"[%d %d %d %d]",N,n1,n2,n);
		//fprintf(stdout,"1.0");
		//exit(0);
		return 0.0;
    } else if (nn <= 1) {
		return 0.0;
    } else if (nn <= HYPERGEO_CACHE_SIZE) {
		return a[nn] ? a[nn] : (a[nn]=gammln(nn+1.0));
    } else {
		return gammln(nn+1.0);
    }
}

double hypergeoD(int N, int n1, int n2, int n) {
    double P=0,c1=0,c2=0,c3=0;
    int i,minIndex;
    minIndex = (n1>n2)?n2:n1;
    c3 = bicolnD(N,n2);
    for (i=n;i<=minIndex;i++) { 
        c1 = bicolnD(n1,i);
        c2 = bicolnD(N-n1,n2-i);
        P += exp((c2+c1)-c3);
    }
	return P;
}

double factlnD(int nn) {
    static double a[HYPERGEO_CACHE_SIZE+1];
    if (nn<0) { 
		fprintf(stderr,"negative factorial\n");
		//fprintf(stdout,"[%d %d %d %d]",N,n1,n2,n);
		fprintf(stdout,"1.0");
		//exit(0);
		return 0.0;
    } else if (nn <= 1) {
		return 0.0;
    } else if (nn <= HYPERGEO_CACHE_SIZE) {
		return a[nn] ? a[nn] : (a[nn]=gammlnD((double)(nn+1.0)));
    } else {
		return gammlnD((double)(nn+1.0));
    }
}
double gammlnD(double xx) {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5396239384953e-5};
    unsigned int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double correlation(double *a, double *b, int n) {
	double x1=0.0,x2=0.0,y1=0.0,y2=0.0,xy=0.0;
	for (int i=0;i<n;i++) {
		x1+=a[i];
		y1+=b[i];
		x2+=a[i]*a[i];
		y2+=b[i]*b[i];
		xy+=a[i]*b[i];
	}
	double numerator = n*xy-x1*y1;
	double denominator = (n*x2 - x1*x1) * (n*y2-y1*y1);
	if (denominator <=0) return -2.0;
	denominator = sqrt(denominator);
	return numerator / denominator;
}

float correlation(float *a, float *b, int n) {
	float x1=0.0,x2=0.0,y1=0.0,y2=0.0,xy=0.0;
	for (int i=0;i<n;i++) {
		x1+=a[i];
		y1+=b[i];
		x2+=a[i]*a[i];
		y2+=b[i]*b[i];
		xy+=a[i]*b[i];
	}
	float numerator = n*xy-x1*y1;
	float denominator = (n*x2 - x1*x1) * (n*y2-y1*y1);
	if (denominator <=0) return -2.0;
	denominator = sqrt(denominator);
	return numerator / denominator;
}

float gammp(float a, float x) {
	float gamser=0.0,gammcf=0.0,gln=0.0;
	if (x < 0.0 || a <= 0.0) {
		return 0;
	}
	if (x < (a+1.0)) {
		gser( a, x, gamser, gln);
		return gamser;
	} else {
		gcf(a, x, gammcf, gln);
		return 1.0-gammcf;
	}
}
float gammq(float a, float x) {
	float gamser = 0.0;
	float gln = 0.0;
	if (x< 0 || a <= 0) {
		fprintf(stderr, "Bad inputs (%f,%f) into gammq\n", a, x);
		return 0.0;
	}
	if (x < a+1.0) {
		gser(a,x,gamser,gln);
		return 1.0-gamser;
	} else {
		gcf(a,x,gamser,gln);
		return gamser;
	}
}

void gser(float a, float x, float &gamser, float &gln) {
	gln = gammln(a);
	gamser =0;
	int itmax = 100;
	float eps = 3.0e-10;
	if (x <= 0.0) {
		gamser = 0;
		return;
	} else {
		float ap = a;
		float del = 1.0/a;
		float sum = del;
		for (int i=1;i<itmax;i++) {
			ap++;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*eps) {
				gamser = sum*exp(-1*x+a*log(x)-gln);
				return;
			}
		}
		fprintf(stderr, "reached max iterations in gser\n");
		return;
	}
}

void gcf(float a, float x, float &gammcf, float& gln) {
	int itmax = 100;
	float eps = 2.0e-8;
	gammcf = 0;
	gln = gammln(a);
	float b = x+1.0-a;
	float fpmin = 1e-30;
	float c = 1.0/fpmin;
	float d = 1.0/b;
	float h = d;
	int i=0;
	float del;
	for (i=1;i<=itmax;i++) {
		float an = i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < fpmin) {
			d = fpmin;
		}
		c = b+an/c;
		if (fabs(c) < fpmin) {
			c = fpmin;
		}
		d = 1.0/d;
		del = d*c;
		h*=del;
		if (fabs(del-1.0) < eps) {
			break;
		}
	}
	if (i>itmax) {
		//fprintf(stderr, "reaced maxit in gcf\n");
	}
	gammcf = exp(-1*x+a*log(x)-gln)*h;
}

float chi2pvalue(float chi2, int df) {
	return gammq(((float)df)/2.0,chi2/2.0);
}
float chi2pvaluelog(float chi2, int df) {
	//return gammqlog(((float)df)/2.0,chi2/2.0);
	return 0;
}



