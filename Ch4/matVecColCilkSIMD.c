/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.3 from
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
 * Implementation of column-wise SIMD nxn matrix-vector multiplication
 * where n is a power of two
 * using Cilk Plus
 * If SIMD not defined, uses array notation
 * If SIMD defined (compile with -DSIMD), uses SIMD pragma
 */
#include <stdio.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <time.h>
#include <float.h>
#include <math.h>

float **allocate2D(int n);
void matvecSerial(float **A, float *x, float *b, int n);

int main(int argc, char **argv){
	if(argc <2){
		fprintf(stderr,"usage: %s exponent\n", argv[0]);
		return 1;
	}
	int exp = strtol(argv[1], NULL, 10);
	int n = 1<<exp;
	float *x = malloc(n*sizeof(float));
	float *b = malloc(n*sizeof(float));
	float *bs = malloc(n*sizeof(float));
	float *temp = malloc(n*sizeof(float));
	float **A = allocate2D(n);
	if(!x || !b || !bs || !temp || !A){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}

	for(int i=0; i<n; i++){
		x[i] = rand();
			for(int j=0; j<n; j++){
				A[i][j] = rand();
			}
	}
	matvecSerial(A, x, bs, n);

	struct timespec tstart,tend;
	float time;
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	#ifdef SIMD
	#pragma SIMD
	for(int i=0; i<n; i++)
		b[i] = 0.0;
	#else
	b[0:n] = 0.0;
	#endif
	for(int i=0;i<n;i++){
	#ifdef SIMD
		#pragma SIMD
		for(int j=0; j<n; j++)
			temp[j] = A[i][j]*x[j];
		for (int k=exp-1; k>=0; k--){
			int k2 = 1<<k;
			#pragma SIMD
			for(int j=0; j<k2; j++)
				temp[j] += temp[j+k2];
		}
	#else
		temp[0:n]=A[i][0:n]*x[0:n];
		for (int k=exp-1; k>=0; k--){
			int k2 = 1<<k;
			temp[0:k2] = temp[0:k2] + temp[k2:k2];
		}
	#endif
		b[i] = temp[0];
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	time = (tend.tv_sec-tstart.tv_sec) + (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("time in s: %f\n", time);
	printf("machine epsilon = %g\n", FLT_EPSILON);
	float max = 0.0;
	for(int j=0;j<n;j++){
		float err = fabs(b[j]-bs[j])/bs[j];
		if(err > max)
			max = err;
	}
	printf("maximum relative difference: %g\n", max);
  return 0;
}

void matvecSerial(float **A, float *x, float *b, int n){
	for(int i=0; i<n; i++)
		b[i] = 0;
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			b[i] += A[i][j]* x[j];
}

float **allocate2D(int n){
		float *temp = malloc(n*n*sizeof(float));
		float **a = malloc(n*sizeof(float *));
		if(!temp || !a)
			return NULL;
		for(int i=0; i<n; i++)
			a[i] = temp + i*n;
		return a;
}
