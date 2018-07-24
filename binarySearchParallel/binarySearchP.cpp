/**
 * Magister en Informatica
 * Universidad Austral de Chile
 * Computacion de Alto Rendimiento
 * INFO335
 * Prof. Dr. Cristobal Navarro
 * 
 * Código base extraído de lab09 de Dr. Héctor Ferrada para Búsqueda 
 * binaria y aplicación de técnicas de HPC.
 * 
 * "Performance comparison of an implementation of a 
 * variant of the binary search algorithm with gap-encoding
 * and search executions parallelization"
 * ------------
 * Alan Keith Paz
 * Due Date: Jul 22, 2018
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <random>
#include <bits/random.h>
#include <time.h>       
#include "include/BasicCDS.h"
#include <omp.h>

using namespace std;
using namespace cds;

#define PRINT 0
#define TEST 1
#define CELLS 64	// minimum number of cells to do bs

ulong REPET = 100000;
uint BLK = 32;	// length of the block for sample array (normal distribution) or for increments (uniform distribution)
uint bM; 		// bits for MAX
uint bMG;		// bits for MAX Gap
uint bMS;		// bits for MAX Sample value
uint p_threads;	// number of threads

/** Structure with all globals parameters program */
typedef struct {
	char prefixResult[300];	// Prefix name of results files in the experiments
	ulong *A;		// original input array
	ulong *X;		// small array for A
	ulong *G;		// gap array for A
	ulong *S;		// sampling array for A and XG
	ulong nWX;		// number of Words for X[]
	ulong nWG;		// number of words for G[]
	ulong nWS;		// number of words for S[]

	ulong *PATT;
	ulong n;
	ulong sizeA, sizeX, sizeG, sizeS;
	ulong min, MAX, MAXG, maxS;

	bool NORMAL;		// flag probability: 1 = NORMAL, 0 = UNIFORM
	uint sigma;		// for normal distribution probability function

	uint s;
	uint sg;		// sampling distance for A[]
	uint bs;
	uint nc;		// number of cells in COUNT
	uint ns;		// number of maximus position of a sample in S
	uint *COUNT;	// array to store the numbers of items per block
	ulong *POS;		// array to store the position in X where each segment begins
} ParProg;

void genArrays(ParProg *par);
void runBS(ParProg *par);
void runAutoBS(ParProg *par);
void runBSSScn(ParProg *par);
void runBSGScn(ParProg *par);
bool binarySearch(ParProg *par, ulong x, ulong *idx);
bool scanBSearch(ParProg *par, ulong x, ulong *idx);
bool gapBSearch(ParProg *par, ulong x, ulong *idx);
void testSearches(ParProg *par);

int main(int argc, char** argv){
	ParProg *par = new ParProg();
	
	/*
	 * Parallel Section
	 * @argv[6] number of threads
	 */
	
	p_threads = atoi(argv[6]);
	omp_set_num_threads(p_threads);

	if(argc < 5){
		cout << "Execution Error! call: ./bsearch <n> <prefixResult> <BLK> <NORMAL flag> [<sigma>] <threads>" << endl;
		exit(EXIT_FAILURE);
	}
	par->n = atoi(argv[1]);
	strcpy(par->prefixResult, "");
	strcpy(par->prefixResult, argv[2]);
	BLK = atoi(argv[3]);
	par->NORMAL = atoi(argv[4]);
	if (par->NORMAL)
		par->sigma = atoi(argv[5]);

	cout << "Parameters..." << endl;
	cout << "n = " << par->n << endl;
	cout << "prefixResult = " << par->prefixResult << endl;
	cout << "NORMAL flag = " << par->NORMAL << endl;
	cout << "Number of Threads = " << p_threads << endl;
	if (par->NORMAL)
		cout << "sigma = " << par->sigma << endl;

	genArrays(par);
	if (TEST){
		testSearches(par);
		cout << " Test OK !! " << endl;
	}

	runAutoBS(par);
	runBS(par);
	runBSSScn(par);
	runBSGScn(par);

	cout << "######################################" << endl;
	return 0;
}

// we look for the final position of val in X[1..n]
void heapify(ulong *X, ulong n, ulong val){
	ulong m=2, i=1;

	while(m<n){
		if (m+1 < n && X[m] < X[m+1])
			m++;

		if(val < X[m]){
			X[i] = X[m];
			i = m;
			m <<= 1;
		}else
			break;
	}
	X[i] = val;
}

// heap sort to sort the array
void sortArray(ulong *X, ulong n){
	ulong i,j,k,val;

	// 1.- create the max-heap...
	for(i=2; i<=n; i++){
		val=X[i];
		k=i/2;
		j=i;
		while(k && val>X[k]){
			X[j]=X[k];
			j=k;
			k=k/2;
		}
		X[j]=val;
	}

	if (TEST){
		for (i=n; i>1; i--){
			if (X[i] > X[i/2]){
				cout << "ERROR. X["<<i<<"]=" << X[i] <<" > X[" << i/2 <<"]=" << X[i/2] << endl;
				exit(1);
			}
		}
	}

	// 2.- Create the final array...
	k=X[n];
	for(i=n-1; i; i--){
		val=X[1];
		heapify(X,i,k);
		k=X[i];
		X[i]=val;
	}

	if (X[1]>X[2]){
		j=X[1];
		X[1] = X[2];
		X[2] = j;
	}

	if(k<X[1])
		X[0]=k;
	else{
		X[0]=X[1];
		X[1]=k;
	}

	if (TEST){
		for (i=1; i<n; i++){
			if (X[i]<X[i-1]){
				cout << "ERROR. X["<<i<<"]=" << X[i] <<" < X[" << i-1 <<"]=" << X[i-1] << endl;
				exit(1);
			}
		}
	}

}

/** 
 * Function that generate X nondecreasing array, PATT array for experiments, sample array for bSearch,
 * G array for gaps and S array for sampling values for gaps
 * @ParProg *par
 */
void genArrays(ParProg *par){
	ulong i, j, k;
	long int num;
	double t, tf;
	char aFile[400];
	char str[100];
	
	par->sizeA = par->n*sizeof(ulong);	// Original size

	if (par->NORMAL){
		par->A = new ulong[par->n+1];

		default_random_engine generator;
		normal_distribution<double> distribution(4*par->sigma, par->sigma);	// (esperanza, varianza)

	    num = distribution(generator);
	    par->A[0] = par->MAX = num;
	// REVISAR
	    for (i=1; i<par->n; i++){
	    	num = distribution(generator);
	    	while (num<0)
	    		num = distribution(generator);
	    	par->A[i] = num;
			if (num > (long int)par->MAX)
				par->MAX = num;
	    }
	    par->A[par->n] = par->A[0];
	    sortArray(par->A, par->n);
	}else{
		par->A = new ulong[par->n];
		par->A[0] = rand()%BLK;
		// REVISAR 
		for (i=1; i<par->n; i++)
			par->A[i] = par->A[i-1] + rand()%BLK;
		par->MAX = par->A[par->n-1];
	}

	par->min = par->A[0];
	cout << "min = " << par->min << ", MAX = " << par->MAX << endl;

	// patterns for experiments
	par->PATT = new ulong[REPET];
	k = par->MAX + BLK;
	// REVISAR SI PODEMOS Y SIRVE PARALELIZAR
	for (i=0; i<REPET; i++)
		par->PATT[i] = rand()%k;
		
	/** Comenzamos a considerar el tiempo de armado aquí */
	t = omp_get_wtime();
	// create sample scan structure...
	par->bs = log2(par->MAX/BLK);
	par->s = pow(2,par->bs);
	par->nc = 1 + par->MAX/par->s;
	par->COUNT = new uint[par->nc];
	par->POS = new ulong[par->nc];
	par->POS[0] = 0;
	uint c=0;
	
	for (i=j=0; i<par->nc; i++){ // REVISAR SI SIRVE Y PODEMOS PARALELIZAR?
		k=(i+1)*par->s;
		for (c=0; j<par->n && par->A[j]<k; j++)
			c++;
		par->COUNT[i]=c;
		if ((i+1)<par->nc)
			par->POS[i+1] = par->POS[i]+c;
	}
	
	/**************************************************************************
	 * create X[] array... 
	 */
	bM = 1+log2(par->MAX);
	par->nWX = (par->n*bM)/(sizeof(ulong)*8);
	if ((par->n*bM)%(sizeof(ulong)*8)>0)
		par->nWX++;

	par->X = new ulong[par->nWX];
	par->sizeX = par->nWX*sizeof(ulong);

	#pragma omp parallel
	{
		#pragma omp for private(i)
		for (i=0; i<par->n; i++)
			setNum64(par->X, i*bM, bM, par->A[i]);
	}

	/**************************************************************
	 * Create array G[] and S[] for gap-encoding
	 * First we search the maximus gap,
	 * then create G[] array for gaps 
	 * and finally create reduce S[] array for sampling
	 **************************************************************/
	ulong mps, gap, maxG;
	maxG = par->A[0];
	#pragma omp parallel
	{
		#pragma omp for private(i,gap)
		for(i=1; i < par->n; i++){
			gap = par->A[i] - par->A[i-1];
			if(gap > maxG)
				maxG = gap;
		}
		#pragma omp barrier // comentar
	}
	par->MAXG = maxG;
	bMG = 1+log2(par->MAXG); // bits 
	par->nWG = (par->n*bMG)/(sizeof(ulong)*8); 
	if ((par->n*bMG)%(sizeof(ulong)*8)>0)
		par->nWG++;

	// create G[]
	par->G = new ulong[par->nWG];
	par->sizeG = par->nWG*sizeof(ulong);
	setNum64(par->G, 0, bMG, par->A[0]);
	#pragma omp parallel
	{
		#pragma omp for private(j,gap)
		for(j=1; j < par->n; j++){
			gap = par->A[j] - par->A[j-1];
			setNum64(par->G, j*bMG, bMG, gap);
		}
	}
	// Calculate value of sg
	par->sg = log2(par->n)*log2(par->n)/16;
	mps = par->n/par->sg;
	par->ns = (uint)mps;
	par->maxS = par->A[par->ns*par->sg];
	bMS = 1+log2(par->maxS);
	par->nWS = (par->ns*bMS)/(sizeof(ulong)*8);
	if ((par->ns*bMS)%(sizeof(ulong)*8)>0)
		par->nWS++;
	// Create S[]
	par->S = new ulong[par->nWS];
	par->sizeS = par->nWS*sizeof(ulong);
	#pragma omp parallel
	{
		#pragma omp for private(k)
		for(k=0; k <= par->ns ; k++){
			setNum64(par->S, k*bMS, bMS, par->A[k*par->sg]);
		}
		#pragma omp barrier // comentar
	}

	tf = omp_get_wtime();
	
	/* Show Create Statistics */
	// For Sample Scan
	cout << "s = " << par->s << endl;
	cout << "bs = " << par->bs << endl;
	cout << "AVG COUNT = " << (float)par->n/(float)par->nc << endl;
	uint aux = par->nc*(sizeof(uint)+sizeof(ulong));
	cout << " Extra size for COUNT[] + POS[] = " << aux/1024.0 << " KiB" << endl;
	// For Reduce X
	cout << "MAX = " << par->MAX << endl;
	cout << "bM = " << bM << endl;
	cout << " size for A[] = " << par->sizeA/(1024.0*1024.0) << " MiB" << endl;
	cout << " size for X[] = " << par->sizeX/(1024.0*1024.0) << " MiB" << endl;
	// For Gap y S
	cout << " size for G[] = " << par->sizeG/(1024.0*1024.0) << " MiB" << endl;
	cout << " size for S[] = " << par->sizeS/(1024.0*1024.0) << " MiB" << endl;
	cout << "Structures Create Time = " << (tf - t) << " Seconds" << endl;
	
	strcpy(aFile, par->prefixResult);
	strcpy(str, "");
	sprintf(str, "CreateTime%d", BLK);
	strcat(aFile, str);
	FILE *fp = fopen(aFile, "a+" );
	// [n] [create time] [threads]
	fprintf(fp, "%ld %lf %d\n", par->n, tf-t, p_threads);
	fclose(fp);
	
	// PRINT DATA STRUCTURES CREATED (FOR TESTING)
	if (PRINT){
		cout << "A[] = ";
		for (i=0; i<par->n; i++)
			cout << par->A[i] << " ";
		cout << endl;

		cout << "X[] = ";
		for (i=j=0; i<par->n; i++, j+=bM)
			cout << getNum64(par->X, j, bM) << " ";
		cout << endl;

		cout << "G[] = ";
		for (i=j=0; i<par->n; i++, j+=bMG)
			cout << getNum64(par->G, j, bMG) << " ";
		cout << endl;

		cout << "S[] = ";
		for (i=j=0; i<(par->n/par->sg); i++, j+=bMS)
			cout << getNum64(par->S, j, bMS) << " ";
		cout << endl;

		cout << "COUNT[0.." << par->nc << "] = ";
		for (i=0; i<par->nc; i++)
			cout << par->COUNT[i] << " ";
		cout << endl;

		cout << "  POS[0.." << par->nc << "] = ";
		for (i=0; i<par->nc; i++)
			cout << par->POS[i] << " ";
		cout << endl;

		c=0;
		for(i=0; c<par->nc; i+=par->COUNT[c], c++);
		if (i!=par->n){
			cout << "ERROR, count cells = " << i << " != n = " << par->n << endl;
			exit(0);
		}
	}
}

int compfunc(const void *a, const void *b){
	return (*(int*)a-*(int*)b);
}

/** 
 * Function that runs the function bsearch() from C
 * @ParProg *par
 */
void runAutoBS(ParProg *par){
	ulong k, nOcc;
	float avgTime;
	double t, tf;
	int *item;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Auto Binary Search on X[] " << endl;

	t = omp_get_wtime();
	#pragma omp parallel
	{
		#pragma omp for private(k,item) reduction(+:nOcc)
		for (k=nOcc=0; k<REPET; k++){
			item = (int*)bsearch(&par->PATT[k], par->A, par->n, sizeof(ulong), compfunc);
			if(item != NULL)
				nOcc += 1;
		}
	}
	tf = omp_get_wtime();
	avgTime = tf - t;
	cout << "Average CPU time per execution: " << (avgTime*1000000.0)/REPET << " Microseconds" << endl;
	cout << "nOcc = " << nOcc << endl;
}

/** paralelizado
 * Function that runs Binary Search
 * @ParProg *par
 */
void runBS(ParProg *par){
	ulong k, nOcc, pos;
	float avgTime;
	char aFile[400];
	char str[100];
	double t, tf;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Binary Search on X[] " << endl;

	t = omp_get_wtime();
	#pragma omp parallel
	{
		#pragma omp for private(k) reduction(+:nOcc)
		for (k=nOcc=0; k<REPET; k++)
			nOcc += binarySearch(par, par->PATT[k], &pos);
	}
	tf = omp_get_wtime();
	avgTime = tf - t;
	cout << "Average CPU time per execution: " << (avgTime*1000000.0)/REPET << " Microseconds" << endl;
	cout << "nOcc = " << nOcc << endl;

	strcpy(aFile, par->prefixResult);
	strcpy(str, "");
	sprintf(str, "bSearchBLK%d", BLK);
	strcat(aFile, str);
	cout << "Resume File: " << aFile << endl;
	// /home/hferrada/Dropbox/UACh/teaching/magister/practicas/busqueda
	FILE *fp = fopen(aFile, "a+" );
	if (par->NORMAL){
		// [n] [REPET] [nOcc] [avg bs-time/exec] [esperanza] [varianza]
		fprintf(fp, "%ld %ld %ld %f %ld %d\n", par->n, REPET, nOcc, (avgTime*1000000.0)/REPET, par->n/2, par->sigma);
	}else{
		// [n] [REPET] [nOcc] [avg bs-time/exec]
		fprintf(fp, "%ld %ld %ld %f\n", par->n, REPET, nOcc, (avgTime*1000000.0)/REPET);
	}
	fclose(fp);
}
/** paralelizado
 * Function that runs Binary Search with Sampling Scan
 * @ParProg *par
 */
void runBSSScn(ParProg *par){
	ulong k, nOcc, pos;
	float avgTime;
	char aFile[400];
	char str[100];
	double t, tf;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Binary Search Sample Scan on X[] " << endl;

	t = omp_get_wtime();
	#pragma omp parallel
	{
		#pragma omp for private(k) reduction(+:nOcc)
		for (k=nOcc=0; k<REPET; k++)
			nOcc += scanBSearch(par, par->PATT[k], &pos);
	}
	tf = omp_get_wtime();
	avgTime = tf - t;
	
	cout << "Average CPU time per execution: " << (avgTime*1000000.0)/REPET << " Microseconds" << endl;
	cout << "nOcc = " << nOcc << endl;

	strcpy(aFile, par->prefixResult);
	strcpy(str, "");
	sprintf(str, "bSearchScanBLK%d", BLK);
	strcat(aFile, str);
	cout << "Resume File: " << aFile << endl;
	FILE *fp = fopen(aFile, "a+" );
	if (par->NORMAL){
		// [n] [REPET] [nOcc] [avg bs-time/exec] [esperanza] [varianza]
		fprintf(fp, "%ld %ld %ld %f %ld %d\n", par->n, REPET, nOcc, (avgTime*1000000.0)/REPET, par->n/2, par->sigma);
	}else{
		// [n] [REPET] [nOcc] [avg bs-time/exec]
		fprintf(fp, "%ld %ld %ld %f\n", par->n, REPET, nOcc, (avgTime*1000000.0)/REPET);
	}
	fclose(fp);
}

/** paralelizado
 * Function that runs Binary Search for Gap encoding with Scan
 * @ParProg *par
 */
void runBSGScn(ParProg *par){
	ulong k, nOcc, pos;
	float avgTime;
	char aFile[400];
	char str[100];
	double t, tf;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Binary Search with Gap-Encoding & Scan on X[] " << endl;
	
	t = omp_get_wtime();
	#pragma omp parallel
	{
		#pragma omp for private(k) reduction(+:nOcc)
		for (k=nOcc=0; k<REPET; k++)
			nOcc += gapBSearch(par, par->PATT[k], &pos);
	}
	tf = omp_get_wtime();
	avgTime = tf - t;
	
	cout << "Average CPU time per execution: " << (avgTime*1000000.0)/REPET << " Microseconds" << endl;
	cout << "nOcc = " << nOcc << endl;

	strcpy(aFile, par->prefixResult);
	strcpy(str, "");
	sprintf(str, "bSearchGapBLK%d", BLK);
	strcat(aFile, str);
	cout << "Resume File: " << aFile << endl;
	FILE *fp = fopen(aFile, "a+" );
	if (par->NORMAL){
		// [n] [REPET] [nOcc] [avg bs-time/exec] [esperanza] [varianza]
		fprintf(fp, "%ld %ld %ld %f %ld %d\n", par->n, REPET, nOcc, (avgTime*1000000.0)/REPET, par->n/2, par->sigma);
	}else{
		// [n] [REPET] [nOcc] [avg bs-time/exec]
		fprintf(fp, "%ld %ld %ld %f\n", par->n, REPET, nOcc, (avgTime*1000000.0)/REPET);
	}
	fclose(fp);
}

/**
 * Function that implements binary search for x on X[]
 * @ParProg * par, @ulong x, @ulong *idx
 * @return 1 for x finded, 0 if not.
 */
bool binarySearch(ParProg *par, ulong x, ulong *idx){
	ulong *A = par->A;
	if (x < A[0] || x >par->MAX)
		return 0;

	ulong l, r, m;

	l=0;
	r=par->n-1;
	m = r/2;

	while (l<=r){
		if (x==A[m]){
			*idx = m;
			return 1;
		}

		if (x<A[m])
			r=m-1;
		else
			l=m+1;

		m=l+(r-l)/2;
	}
	return 0;
}

/**
 * Function that implements binary search with scan for x on X[]
 * @ParProg * par, @ulong x, @ulong *idx
 * @return 1 for x finded, 0 if not.
 */
bool scanBSearch(ParProg *par, ulong x, ulong *idx){
	ulong *X = par->X;
	if (x < par->min || x >par->MAX)
		return 0;

	uint pos = x>>par->bs;
	ulong c=par->COUNT[pos];
	ulong m, xm, l=par->POS[pos];

	if (c > CELLS){
		// Binary searching in the X-segment...
		ulong r;

		r=l+c-1;
		m=(l+r)>>1;
		xm = getNum64(X, m*bM, bM);

		while (l<=r){
			if (x<xm)
				r=m-1;
			else{
				if (x==xm){
					*idx = m;
					return 1;
				}
				l=m+1;
			}
			m=(l+r)>>1;
			xm=getNum64(X, m*bM, bM);
		}
	}else{
		// Scanning a X-segment of maximum c cells...
		// for (m=l*bM, xm=getNum64(X, m, bM); xm<x; m+=bM, xm=getNum64(X, m, bM))
		m=l*bM;
		xm=getNum64(X, m, bM);
		while (xm<x){
			m+=bM;
			xm=getNum64(X, m, bM);
		}

		if (xm==x){
			*idx = m/bM;
			return 1;
		}
	}

	return 0;
}

/**
 * Function that implements binary search for gap-encoding for x on X[]
 * @ParProg * par, @ulong x, @ulong *idx
 * @return 1 for x finded, 0 if not.
 */
bool gapBSearch(ParProg *par, ulong x, ulong *idx){
	ulong *G = par->G;
	ulong *S = par->S;
	ulong l, r, m, xm;
	if (x < par->min || x >par->MAX)
		return 0;

	l = 0;
	r = par->ns;
	m = r>>1;
	xm = getNum64(S, m*bMS, bMS);
	// We perform a binary search in S looking for the segment where x is
	if(x<= par->maxS){
		while (l<=r){
			if(x<xm)
				r=m-1;
			else{
				if(x==xm){
					*idx = m*par->sg;
					return 1;
				}
				l=m+1;
			}
			m=(l+r)>>1;
			xm=getNum64(S, m*bMS, bMS);
		}
	}
	else
		m = par->ns;
	l = m*par->sg;
	xm = getNum64(S, m*bMS, bMS);

	m = l*bMG;
	while (xm<x){
		m+=bMG;
		xm += getNum64(G, m, bMG);
	}
	if (xm==x){
		*idx = m/bMG;
		return 1;
	}// */
	return 0;
}

/**
 * Function that test the integrity of results from the 3 differents search methods
 * @ParProg * par
 */
void testSearches(ParProg *par){
	ulong k, p1=par->n+1, p2=par->n+2, p3=par->n+3;
	bool found;

	for (k=0; k<REPET; k++){
		p1=par->n+1;
		p2=0;
		p3=0;

		found = binarySearch(par, par->PATT[k], &p1);
		if (scanBSearch(par, par->PATT[k], &p2))
			found = 1;
		if (gapBSearch(par, par->PATT[k], &p3))
			found = 1;

		if (found && (par->A[p1] != par->A[p2] || par->A[p1] != par->A[p3])){
			cout << "ERROR. patt[" <<k<< "] = " << par->PATT[k] << " binarySearch().pos = " << p1 << " != scanBSearch().pos = " << p2 << " != gapBSearch().pos = " << p3 << endl;
			exit(1);
		}
	}
}
