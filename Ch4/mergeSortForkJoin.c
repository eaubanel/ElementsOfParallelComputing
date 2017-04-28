/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.6 from
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
 * Fork-join merge sort, using Cilk Plus. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

/*sort elements of a with index i in [lower .. upper) into array b
 * in parallel using Cilk Plus
 * array a is modified
 */
void parMergeSort(int *a, int lower, int upper, int *b);
//parallel merge using Cilk Plus of sorted subarrays
// i in [low1..up1) and i in [low2..up2) at index start
void parMerge(int *a, int low1, int up1, int low2, int up2, int *b, int start);
//first index in [low..up) such that a[index] > a[ikey]
int binarySearch(int *a, int low, int up, int ikey);
//merge a[low1 .. up1) and a[low2 .. up2) into b at index start
void sequentialMerge(int *a, int low1, int up1, int low2, int up2, int *b, int start);
//comparison function for sequential sort
int comparefunc(const void * arg1, const void * arg2);
//sequential merge sort
void mergeSort(int *a, int lower, int upper, int *b);
//merge a[lower .. mid) and a[mid .. upper) into b
void merge(int *a, int lower, int mid, int upper, int *b);

int cutoff; //cutoff for recursion

int main(int argc, char **argv){
	int *a; //array to be sorted
	int *b; //array used for merging
	int *bs; //array used for sequential sort
	int n; //size of arrays

	struct timespec tstart,tend; 
  float timer;

	if(argc <3){
		fprintf(stderr,"usage: %s n cutoff\n", argv[0]);
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	cutoff = strtol(argv[2], NULL, 10);
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
	parMergeSort(a, 0, n, b);
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

void parMergeSort(int *a, int lower, int upper, int *b){
	if(upper - lower < cutoff)
		mergeSort(a, lower, upper, b);
	else {
		int mid = (upper + lower)/2;
		cilk_spawn parMergeSort(b, lower, mid, a);
		parMergeSort(b, mid, upper, a);
		cilk_sync;
		parMerge(a, lower, mid, mid, upper, b, lower);
	}
}

void parMerge(int *a, int low1, int up1, int low2, int up2, int *b, int start){
	int k1 = up1 - low1;
	int k2 = up2 - low2;
	int mid1, mid2;
	if(k1+k2 < cutoff)
		sequentialMerge(a, low1, up1, low2, up2, b, start);
	else{
		if(k1 >= k2){
			mid1 = (low1+up1-1)/2;
			mid2 = binarySearch(a, low2, up2, mid1);
		} else{
			mid2 = (low2+up2-1)/2;
			mid1 = binarySearch(a, low1, up1, mid2)-1;
			mid2++;
		}
		cilk_spawn parMerge(a, low1, mid1+1, low2, mid2, b, start);
		start = start + mid1 - low1 + 1 + mid2 - low2;
		parMerge(a, mid1+1, up1, mid2, up2, b, start);
		cilk_sync;
	}
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
