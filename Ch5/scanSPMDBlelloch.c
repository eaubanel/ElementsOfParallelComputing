/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 5.1
 * (with Blelloch scan) from
 * Elements of Parallel Computing, by Eric Aubanel, 2016, CRC Press.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------
 * Implementation of SPMD inclusive prefix sum of n integers and p threads 
 * using OpenMP, where n divisible by p, and Blelloch scan used 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

//Performs inclusive prefix sum of array a of length n
void prefixSum(int *a, int n);
//Performs parallel (OpenMP) exclusive prefix sum
//of p values (one per thread, where power of 2 threads) 
//using Blelloch algorithm. Returns idth element.
//Contains orphaned OpenMP directives (barrier)
int parPrefixSumBlelloch(int *a, int p, int id);
int isPowerOf2(int n);
	
int main(int argc, char **argv){
	if(argc < 2){
		fprintf(stderr,"usage: %s n\n", argv[0]);
		return 1;
	}
	int n = strtol(argv[1], NULL, 10);
	int p = omp_get_max_threads();
	if(!isPowerOf2(p)){
		fprintf(stderr,"p (%d) must be power of 2\n", p);
		return 1;
	}
	if(n%p){
		fprintf(stderr,"n must be divisible by p (%d)\n", p);
		return 1;
	}
	int np = n/p; //number of elements per thread
	int *a = malloc(n*sizeof(int));
	int *as = malloc(n*sizeof(int)); //for sequential verification
	if(!a || !as){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	for(int i=0; i<n; i++)
		a[i] = rand()%10;
	memcpy(as, a, n*sizeof(int));
	// sequential prefix sum for verification
	prefixSum(as, n);

	struct timespec tstart, tend;
	float time;
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	int b[p];
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int start = id*np;
		prefixSum(a+start, np);	
		b[id] = a[(id+1)*np-1]; //bug in Algorithm 5.1
		#pragma omp barrier
		//modified from Algorithm 5.1, 
		//since exclusive prefix sum of b used
		int t = parPrefixSumBlelloch(b, p, id);
		for(int j=0; j<np; j++)
			a[start+j] += t;
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	time = (tend.tv_sec-tstart.tv_sec) + (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("parallel time in s: %f\n", time);
	int passed = 1;
	for(int i=0; i<n; i++)
		if(a[i] != as[i]){
			fprintf(stderr, "a[%d]=%d, as[%d] = %d\n", i, a[i], i, as[i]);
			passed = 0;
		}
	if(passed)
		printf("result verified\n");
	return 0;
}

void prefixSum(int *a, int n){
	int sum = 0;
	for(int i=0; i<n; i++){
		sum += a[i];
		a[i] = sum;
	}
}

int parPrefixSumBlelloch(int *a, int p, int id){
	int logp = 0;
	//first pass up tree (reduction)
	for(int j=1; j<p; j<<=1, logp++){
		int tj = j<<1;
		if(0 == id%tj) 
			a[id+tj-1] += a[id+j-1];
		#pragma omp barrier
	}
	//only thread 0 in first iteration of next loop.
	if(!id)
		a[p-1] = 0;
	//second pass down tree
	for(int j=1<<(logp-1); j>=1; j>>=1){
		int tj = j<<1;
		if(0 == id%tj){
			int t = a[id+j-1];
			a[id+j-1] = a[id+tj-1];
			a[id+tj-1] = a[id+j-1] + t;
		}	
		#pragma omp barrier
	}
	return a[id];
}

int isPowerOf2(int n){
	while(n){
		if(n & 1)
			break;
		n >>= 1;
	}
	return (1 == n? 1:0);
}
