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
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "include/BasicCDS.h"

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

// Structure with all globals parameters program
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
	uint ns;		// number of elements in S
	uint *COUNT;	// array to store the numbers of items per block
	ulong *POS;		// array to store the position in X where each segment begins
} ParProg;

void genArrays(ParProg *par);
void runAutoBS(ParProg *par);
void runBS(ParProg *par);
void runBSSScn(ParProg *par);
void runBSGScn(ParProg *par);
bool binarySearch(ParProg *par, ulong x, ulong *idx);
bool scanBSearch(ParProg *par, ulong x, ulong *idx);
bool gapBSearch(ParProg *par, ulong x, ulong *idx);
void testSearches(ParProg *par);

// 50 /home/hferrada/Dropbox/UACh/teaching/magister/HPC/slides/practicaBSearch/results/ 8 1 10
// 100000000 /home/hferrada/Dropbox/UACh/teaching/magister/HPC/slides/practicaBSearch/results/ 40000 1 10000000
int main(int argc, char** argv){
	ParProg *par = new ParProg();

	if(argc < 5){
		cout << "Execution Error! call: ./bsearch <n> <prefixResult> <BLK> <NORMAL flag> [<sigma>]" << endl;
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

	cout << "##################" << endl;
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
	/*cout << "H[] =" << endl;
	for (i=1; i<=n; i++)
		cout << X[i] << " ";
	cout << endl;*/

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

// generate X nondecreasing array, PATT array for experiments, sample array for bSearch,
// G array for gaps and S array for sampling values for gaps
void genArrays(ParProg *par){
	ulong i, j, k;
	long int num;
	char aFile[400];
	char str[100];
	clock_t t;
	float totaltime;

	par->sizeA = par->n*sizeof(ulong);	// Original size

	if (par->NORMAL){
		par->A = new ulong[par->n+1];

		default_random_engine generator;
		normal_distribution<double> distribution(4*par->sigma, par->sigma);	// (esperanza, varianza)

	    num = distribution(generator);
	    par->A[0] = par->MAX = num;

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
		for (i=1; i<par->n; i++)
			par->A[i] = par->A[i-1] + rand()%BLK;
		par->MAX = par->A[par->n-1];
	}

	par->min = par->A[0];
	cout << "min = " << par->min << ", MAX = " << par->MAX << endl;

	// patterns for experiments
	par->PATT = new ulong[REPET];
	k = par->MAX + BLK;
	for (i=0; i<REPET; i++)
		par->PATT[i] = rand()%k;

	// create sample scan structure...
	// We take time from here...
	t = clock();
	par->bs = log2(par->MAX/BLK);
	par->s = pow(2,par->bs);
	par->nc = 1 + par->MAX/par->s;
	par->COUNT = new uint[par->nc];
	par->POS = new ulong[par->nc];
	par->POS[0] = 0;
	uint c=0;
	for (i=j=0; i<par->nc; i++){
		k=(i+1)*par->s;
		for (c=0; j<par->n && par->A[j]<k; j++)
			c++;
		par->COUNT[i]=c;
		if ((i+1)<par->nc)
			par->POS[i+1] = par->POS[i]+c;
	}

	// **************************************************************************
	// create X[] array...
	bM = 1+log2(par->MAX);
	par->nWX = (par->n*bM)/(sizeof(ulong)*8);
	if ((par->n*bM)%(sizeof(ulong)*8)>0)
		par->nWX++;

	par->X = new ulong[par->nWX];
	par->sizeX = par->nWX*sizeof(ulong);

	for (i=j=0; i<par->n; i++, j+=bM)
		setNum64(par->X, j, bM, par->A[i]);

	/**************************************************************
	 * Create array G[] and S[] for gap-encoding
	 * First we search the maximus gap,
	 * then create G[] array for gaps 
	 * and finally create reduce S[] array for sampling
	 **************************************************************/
	ulong mps, gap, maxG;
	// max gap
	maxG = par->A[0];
	for(i=1; i < par->n; i++){
		gap = par->A[i] - par->A[i-1];
		if(gap > maxG)
			maxG = gap;
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
	for(j=1; j < par->n; j++){
		gap = par->A[j] - par->A[j-1];
		setNum64(par->G, j*bMG, bMG, gap);
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

	for(k=0; k <= par->ns ; k++){
		setNum64(par->S, k*bMS, bMS, par->A[k*par->sg]);
	}
	t = clock() - t;
	totaltime = (float)t/CLOCKS_PER_SEC;
	/* Mostrar estadísticas de armado */
	// Para Sample Scan
	cout << "s = " << par->s << endl;
	cout << "bs = " << par->bs << endl;
	cout << "AVG COUNT = " << (float)par->n/(float)par->nc << endl;
	uint aux = par->nc*(sizeof(uint)+sizeof(ulong));
	cout << " Extra size for COUNT[] + POS[] = " << aux/1024.0 << " KiB" << endl;
	// Para Reduce X
	cout << "MAX = " << par->MAX << endl;
	cout << "bM = " << bM << endl;
	cout << " size for A[] = " << par->sizeA/(1024.0*1024.0) << " MiB" << endl;
	cout << " size for X[] = " << par->sizeX/(1024.0*1024.0) << " MiB" << endl;
	// Para Gap y S
	cout << " size for G[] = " << par->sizeG/(1024.0*1024.0) << " MiB" << endl;
	cout << " size for S[] = " << par->sizeS/(1024.0*1024.0) << " MiB" << endl;
	cout << "Structures Create Time = " << totaltime << " Seconds" << endl;
	
	strcpy(aFile, par->prefixResult);
	strcpy(str, "");
	sprintf(str, "CreateTime%d", BLK);
	strcat(aFile, str);
	FILE *fp = fopen(aFile, "a+" );
	// [n] [create time]
	fprintf(fp, "%ld %lf\n", par->n, totaltime);
	fclose(fp);

	// PRINT DATA STRUCTURES CREATED
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

		/*cout << "PATT[] =" << endl;
		for (i=0; i<REPET; i++)
			cout << par->PATT[i] << " ";
		cout << endl;*/

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

void runAutoBS(ParProg *par){
	ulong k, nOcc;
	float avgTime;
	clock_t t;
	int *item;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Auto Binary Search on X[] " << endl;

	t = clock();
	for (k=nOcc=0; k<REPET; k++){
		item = (int*)bsearch(&par->PATT[k], par->A, par->n, sizeof(ulong), compfunc); // Esto se puede hacer con X y el tamaño de x
		if(item != NULL)
			nOcc += 1;
	}
	
	t = clock() - t;
	avgTime = (float)t/CLOCKS_PER_SEC;
	cout << "Average CPU time per execution: " << (avgTime*1000000.0)/REPET << " Microseconds" << endl;
	cout << "nOcc = " << nOcc << endl;
}

void runBS(ParProg *par){
	ulong k, nOcc, pos;
	float avgTime;
	char aFile[400];
	char str[100];
	clock_t t;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Binary Search on X[] " << endl;

	t = clock();
	for (k=nOcc=0; k<REPET; k++)
		nOcc += binarySearch(par, par->PATT[k], &pos);

	t = clock() - t;
	avgTime = (float)t/CLOCKS_PER_SEC;
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

void runBSSScn(ParProg *par){
	ulong k, nOcc, pos;
	float avgTime;
	char aFile[400];
	char str[100];
	clock_t t;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Binary Search Sample Scan on X[] " << endl;

	t = clock();
	for (k=nOcc=0; k<REPET; k++)
		nOcc += scanBSearch(par, par->PATT[k], &pos);

	t = clock() - t;
	avgTime = (float)t/CLOCKS_PER_SEC;
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

void runBSGScn(ParProg *par){
	ulong k, nOcc, pos;
	float avgTime;
	char aFile[400];
	char str[100];
	clock_t t;

	cout << "_________________________________________________" << endl;
	cout << "  Executing " << REPET << " Binary Search with Gap-Encoding & Scan on X[] " << endl;
	
	t = clock();

	for (k=nOcc=0; k<REPET; k++)
		nOcc += gapBSearch(par, par->PATT[k], &pos);

	t = clock() - t;
	avgTime = (float)t/CLOCKS_PER_SEC;
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

// binary search for x on X[]
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

bool gapBSearch(ParProg *par, ulong x, ulong *idx){
	ulong *G = par->G;
	ulong *S = par->S;
	if (x < par->min || x >par->MAX)
		return 0;

	ulong l, r, m, xm, aux; //, gm,i;
	/* boolean flagR for right or left shift in G[]'s scan
	 * true if it's right
	 * false if it's left
	 */
	//bool flagR;
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
		aux = r;
	}
	else
		aux= par->ns;
	l = aux*par->sg;
	xm = getNum64(S, aux*bMS, bMS);

	// Scan through G to find if x is in the array or not.
	m = l*bMG;
	while (xm<x){
		m+=bMG;
		xm += getNum64(G, m, bMG);
	}
	if (xm==x){
		*idx = m/bMG;
		return 1;
	}
	
	/* Binary Search for G // No funcional, pues se deben ocupar 
	 * muchos for para poder ir actualizando m.... 
	 * Doble de tiempo en ejecución!!
	//m = (l+r)>>1; 
	//gm = 0;
	//for (i = l+1; i <= m; i++){
	//	gm += getNum64(G, i*bMG, bMG);
	//}
	xm = gm + getNum64(S, aux*bMS, bMS); // calcular m con el sample y las diferencias
	// We perform now a scan on G[] 
	 * while (l<=r){
		cout << "X = " << x << "; XM = " << xm << " --- ";
		if(x<xm){
			r=m-1;
			flagR = true;
		}
		else{
			if(x==xm){
				*idx = m;
				return 1;
			}
			l=m+1;
			flagR = false;
		}
		m=(l+r)>>1;
		gm = 0;
		if(flagR){
			for (i = m+1; i <= r; i++){
				gm -= getNum64(G, i*bMG, bMG);
				cout << "ctm" << i;
			}
		}
		else{
			for (i = l; i <= m; i++){
				gm += getNum6998434(G, i*bMG, bMG);
				cout << "ctm2" << i;
			}
		}
		xm += gm; // actualizar m con sample y las diferencias
		//cout << "sumado";
	}*/
	return 0;
}

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
