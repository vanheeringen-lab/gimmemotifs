#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef STATISTICS_H
#define STATISTICS_H


#define MAXMEMORY 1e7
#define EMPTY_STAT 1000.0

class StatMemory {
public:
	int length;
	int maxBytes;
	int curBytes;
	float denominator;
	int N;
	int n1;
	float** result;
	StatMemory(int maxsize, int NN, int nn1);
	~StatMemory();
	float getStat(int n1, int n);
};

float hypergeo(int N, int n1, int n2, int n);
float loghypergeo(int N, int n1, int n2, int n);
float iloghypergeo(int N, int n1, int n2, int n);
// for logbinomial, n=num trial, k=successes, r= rate of suceess, N= size of sample to determin rate
float logbinomial(int n, int k, float r, int N);
double hypergeoD(int N, int n1, int n2, int n);
float ttest(float* data, int n1, float* data2, int n2, float &t);
float rankSumStat(int ranksum, int n1, int N);
double logPoisson(int x, double lambda);
double cumulativePoisson(int x, double lambda);
double logCumulativePoisson(int x, double lambda);
double ilogCumulativePoisson(int x, double lambda);

double loghypergeoD(unsigned int N, unsigned int n1, unsigned int n2, unsigned int n);
double iloghypergeoD(unsigned int N, unsigned int n1, unsigned int n2, unsigned int n);
// for logbinomial, n=num trial, k=successes, r= rate of suceess, N= size of sample to determin rate
double logbinomialD(unsigned int n, unsigned int k, double r, unsigned int N);
double ilogbinomialD(unsigned int n, unsigned int k, double r, unsigned int N);


float gammln(float xx);
double gammlnD(double xx);
void avevar(float* data, int n, float &ave, float &var);
float betai(float a, float b, float x);
float logbetai(float a, float b, float x);
double logbetai(double a, double b, double x);
float betacf(float a, float b, float x);
double betacfD(double a, double b, double x);
float bicoln(int n, int k);
float factln(int n);
double bicolnD(unsigned int n, unsigned int k);
double factlnD(int n);
float correlation(float* a, float* b, int n);
double correlation(double* a, double* b, int n);
void gcf(float a, float x, float &gammcf, float& gln);
float chi2pvalue(float chi2, int df);
void gser(float a, float x, float &gamser, float &gln);

#define HYPERGEO_CACHE_SIZE 1000000

#endif
