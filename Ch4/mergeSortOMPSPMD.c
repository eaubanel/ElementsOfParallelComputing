/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.13 from
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
 * Merge sort of n ints, using SPMD style OpenMP with nt threads,
 * where n mod nt = 0 and nt = 2^j 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

//parallel merge sort using SPMD style Openmp (nt threads)
//where n mod nt = 0 and nt must be a power of 2
//b is initialized with values of a
//returns pointer to sorted array
int *parMergeSort(int *a, int *b, int n);
//nmt threads merge a[low1..low2) with a[low2..up2] into b at index low1
void spmdMerge(int *a, int low1, int low2, int up2, int *b, int n, int nmt, int id, int nt);
//first index in [low..up) such that a[index] > a[ikey]
int binarySearch(int *a, int low, int up, int ikey);
//merge a[low1 .. up1) and a[low2 .. up2) into b at index start
void sequentialMerge(int *a, int low1, int up1, int low2, int up2, int *b, int start);
//returns 1 if n is a power of 2, 0 otherwise
int isPowerOf2(int n);
//swap two int* pointers
void swap(int **a, int **b);
//comparison function for sequential sort
int comparefunc(const void * arg1, const void * arg2);
//sequential merge sort of a[lower..upper) into b
void mergeSort(int *a, int lower, int upper, int *b);
//merge a[lower .. mid) and a[mid .. upper) into b
void merge(int *a, int lower, int mid, int upper, int *b);


int main(int argc, char **argv){
	int *a; //array to be sorted
	int *b; //array used for merging
	int *bs; //array used for sequential sort
	int n; //size of arrays

	struct timespec tstart,tend; 
  float timer;

	if(argc < 2){
		fprintf(stderr,"usage: %s n\n", argv[0]);
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	int nt = omp_get_max_threads();
	if(!isPowerOf2(nt)){
		fprintf(stderr, "must have power of 2 number of threads\n");
		return 1;
	}
	if(n%nt != 0){
		fprintf(stderr, "%d must be divisible by number of threads\n", n);
		return 1;
	}
	a = malloc(n*sizeof(int));
	b = malloc(n*sizeof(int));
	bs = malloc(n*sizeof(int));
	if(a == NULL || b == NULL || bs == NULL){
	      	fprintf(stderr,"couldn't allocate array of %d ints\n", n);
		return 1;
	}
	srand(time(NULL));
	for(int i=0; i<n; i++)
		a[i] = rand();

	memcpy(bs, a, n*sizeof(int)); //copy array a to bs
	//sequential sort
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	qsort(bs, n, sizeof(int), comparefunc);
	clock_gettime(CLOCK_MONOTONIC, &tend);
  timer = (tend.tv_sec-tstart.tv_sec) +
        (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("sequential time in s: %f\n", timer);

	//parallel sort
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	memcpy(b, a, n*sizeof(int));
	b = parMergeSort(a, b, n);
	clock_gettime(CLOCK_MONOTONIC, &tend);
  timer = (tend.tv_sec-tstart.tv_sec) +
        (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("parallel time in s: %f\n", timer);

	int passed = 1;
	for(int i=0; i<n; i++)
		if(b[i] != bs[i]){
			printf("i=%d: b[i]=%d, bs[i]=%d\n", i, b[i], bs[i]); 
			passed = 0;
		}
	if(passed)
		printf("result verified\n");
  return 0;
}

int *parMergeSort(int *a, int *b, int n){
	int nt = omp_get_max_threads();
	int temp = nt;
	int lognt = 0;
	while(temp >>= 1)
		lognt++;
	#pragma omp parallel firstprivate(a, b) 
	{
		int id = omp_get_thread_num();
		int lower = id*n/nt;
		int upper = (id+1)*n/nt;
		//each thread sorts its chunk
		mergeSort(a, lower, upper, b);
		#pragma omp barrier
		int nmt = 1;
		for(int i=1; i<=lognt; i++){
			swap(&a, &b);
			int chunk = nmt*n/nt;
			nmt *= 2;
			int idc = (id/nmt)*nmt;
			int low1 = idc*n/nt;
			int low2 = low1 + chunk;
			int up2 = low2 + chunk - 1;
			spmdMerge(a, low1, low2, up2, b, n, nmt, id, nt);
			#pragma omp barrier
		}
	}
	if(lognt%2)
		return a;
	else 
		return b;
}

void spmdMerge(int *a, int low1, int low2, int up2, int *b, int n, int nmt, int id, int nt){
	int idm = id%nmt;
	int lowX = idm*n/(2*nt) + low1;
	int upX = (idm+1)*n/(2*nt) + low1 - 1;
	int lowY, upY;
	if(idm != 0)
		lowY = binarySearch(a, low2, up2+1, lowX-1);
	else
		lowY = low2;
	if(idm < nmt-1)
		upY = binarySearch(a, lowY, up2+1, upX) - 1;
	else
		upY = up2;
	int start = lowX + lowY - low2;
	sequentialMerge(a, lowX, upX+1, lowY, upY+1, b, start);
}

int binarySearch(int *a, int low, int up, int ikey){
	up--; //up now refers to index of last element
	int key = a[ikey];
	while(low <= up){
		int mid = (low+up)/2;
		if(a[mid] == key)
			low = mid+1;
		else if (a[mid] > key)
			up = mid-1;
		else
			low = mid+1;
	}
	if(up < low)
		return low;
	else if(a[up] <= key)
		return up+1;
	else 	
		return up;
}

void sequentialMerge(int *a, int low1, int up1, int low2, int up2, int *b, int start){
	int i = low1, j = low2;
	int nel = up1-low1+up2-low2;
	for(int k=start; k<start+nel; k++){
		if((i < up1) && ((j >= up2) || (a[i] <= a[j])))
			b[k] = a[i++];
		else
			b[k] = a[j++];
	}
}

int isPowerOf2(int n){
	while(n){
		if(n & 1)
			break;
		n >>= 1;
	}
	return (1 == n? 1:0);
}

void swap(int **a, int **b){
	int *temp = *a;
	*a = *b;
	*b = temp;
}


int comparefunc(const void * arg1, const void * arg2){
  const int x = * (const int *)arg1; 
  const int y = * (const int *)arg2; 
	return x - y;
}

void mergeSort(int *a, int lower, int upper, int *b){
	if(upper - lower < 2)
		return;
	int mid = (upper + lower)/2;
	mergeSort(b, lower, mid, a);
	mergeSort(b, mid, upper, a);
	merge(a, lower, mid, upper, b);
}

void merge(int *a, int lower, int mid, int upper, int *b){
	int i = lower, j = mid;
	for(int k=lower; k<upper; k++)
		if((i < mid) && ((j >= upper) || (a[i] <= a[j])))
			b[k] = a[i++];
		else
			b[k] = a[j++];
}
