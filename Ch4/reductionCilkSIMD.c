/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.2 from
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
 * Implementation of divergence-free reduction using Cilk Plus
 * For an array of 2^n integers
 * If SIMD not defined, uses array notation
 * If SIMD defined (compile with -DSIMD), uses SIMD pragma
 */
#include <stdio.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <time.h>

int main(int argc, char **argv){
	if(argc <2){
		fprintf(stderr,"usage: %s exponent\n", argv[0]);
		return 1;
	}
	int exp = strtol(argv[1], NULL, 10);
	int n = 1<<exp;
	int *a = malloc(n*sizeof(int));
	if(!a){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	struct timespec tstart, tend;
	float time;

	for(int i=0; i<n; i++)
		a[i] = rand()%10;
	// sequential reduction for verification
	int sum = 0;
	for(int i=0; i<n; i++)
		sum += a[i];

	clock_gettime(CLOCK_MONOTONIC, &tstart);
	for (int k=exp-1; k>=0; k--) {
		int j = 1<<k;
	#ifdef SIMD
		#pragma simd
		for(int i=0; i<j; i++)
			a[i] += a[i+j];	
	#else
		a[0:j] = a[0:j]+a[j:j];
	#endif
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	time = (tend.tv_sec-tstart.tv_sec) + (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("SIMD time in s: %f\n", time);
	if(a[0] == sum)
		printf("result verified\n");
	else
		printf("parallel sum: %d, sequential sum: %d\n", a[0], sum);
	return 0;
}
