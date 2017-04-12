/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithms 3.2 and 3.3 from
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
 * Recursive merge sort. Merges alternatively between 2 arrays. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*sort elements of a with index i in [lower .. upper) into array b
 * array a is modified
 */
void mergeSort(int *a, int lower, int upper, int *b);
//merge a[lower .. mid) and a[mid .. upper) into b
void merge(int *a, int lower, int mid, int upper, int *b);

int main(int argc, char **argv){
	int *a; //array to be sorted
	int *b; //array used for merging
	int n; //size of arrays

	if(argc <2){
		fprintf(stderr,"usage: %s num_points\n", argv[0]);
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	a = malloc(n*sizeof(int));
	b = malloc(n*sizeof(int));
	if(a == NULL || b == NULL){
	      	fprintf(stderr,"couldn't allocate array of %d ints\n", n);
		return 1;
	}
	for(int i=0; i<n; i++)
		if(scanf("%d", &a[i]) != 1){
			fprintf(stderr,"missing/invalid data\n");
			return 1;
		}

	//copy array a to b
	memcpy(b, a, n*sizeof(int));
	mergeSort(a, 0, n, b);

	for(int i=0; i<n; i++)
		printf("%d\n", b[i]);
}

void mergeSort(int *a, int lower, int upper, int *b){
	if(upper - lower < 2)
		return;
	int mid = (upper + lower)/2;
	mergeSort(b, lower, mid, a);
	mergeSort(b, mid, upper, a);
	merge(a, lower, mid, upper, b);
}

//bug in Alg. 3.3, which redundantly initializes and updates k
void merge(int *a, int lower, int mid, int upper, int *b){
	int i = lower, j = mid;
	for(int k=lower; k<upper; k++)
		if((i < mid) && ((j >= upper) || (a[i] <= a[j])))
			b[k] = a[i++];
		else
			b[k] = a[j++];
}
