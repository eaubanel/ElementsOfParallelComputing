/* Copyright 2017 Eric Aubanel
 * This file contains code implementing SIMD subset sum from Chap. 4 of
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
 * Implementation of T/F subset sum problem using dynamic programming
 * using SIMD with Cilk Plus
 * If SIMD not defined, uses array notation
 * If SIMD defined (compile with -DSIMD), uses SIMD pragma
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cilk/cilk.h>

int main(int argc, char **argv){
  int R; // max magnitude of elements in set
  int n; // number of points in set
	int S; // target sum = nR/4
	int *s; // set of points (starts at index 1)
	char *Fs; // Dynamic programming table (sequential)
	int *F; // Dynamic programming table

	struct timespec tstart,tend; 
  float timer;

	if(argc < 3){
		fprintf(stderr,"usage: %s R n [seed]\n", argv[0]);
		return 1;
	}
	R = strtol(argv[1], NULL, 10);
	n = strtol(argv[2], NULL, 10);
  S = n*R/4;
  int m = S+1;
  s = malloc((n+1)*sizeof(int));
  Fs = calloc((n+1)*m,sizeof(char));
  F = calloc((n+1)*m,sizeof(int));
	if(!s || !F || !Fs){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	if(4 == argc) 
		srand(strtol(argv[3], NULL, 10));
	else
		srand(time(NULL));
  for(int i = 1; i<=n;i++){
		s[i] = rand()%R;
		printf("%d ", s[i]);
		F[i*m] = 1;
		Fs[i*m] = 1;
  }
	printf("\n sum = %d\n", S);
  F[m+s[1]] = 1;
  Fs[m+s[1]] = 1;
	//sequential solution
  for(int i=2;i<=n;i++){
    for(int j=1; j<s[i];j++){
      Fs[i*m+j] = Fs[(i-1)*m+j];
    }
    for(int j=s[i];j<=S;j++){
      Fs[i*m+j] = Fs[(i-1)*m+j] || Fs[(i-1)*m+j-s[i]];
    }
  }

	/* SIMD solution
	* change || to + so that second loop can be vectorized
	*/
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	for(int i=2;i<=n;i++){
		#ifdef SIMD
		#pragma SIMD
    for(int j=1; j<s[i];j++){
      F[i*m+j] = F[(i-1)*m+j];
    }
		#pragma SIMD
    for(int j=s[i];j<=S;j++){
      F[i*m+j] = F[(i-1)*m+j] + F[(i-1)*m+j-s[i]];
    }
		#else
		F[i*m+1:s[i]-1] = F[(i-1)*m+1:s[i]-1];
    F[i*m+s[i]:S-s[i]+1] = F[(i-1)*m+s[i]:S-s[i]+1] + F[(i-1)*m:S-s[i]+1];
		#endif
  }
	clock_gettime(CLOCK_MONOTONIC, &tend);
  timer = (tend.tv_sec-tstart.tv_sec) +
        (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
  printf("%s\n",F[n*m+S]?"true":"false");
	printf("time in s: %f\n", timer);

	//verify
	int passed = 1;
	for(int i=2;i<=n;i++){
    for(int j=1; j<s[i];j++){
			if((F[i*m+j] && !Fs[i*m+j]) || (!F[i*m+j] && Fs[i*m+j])){
				printf("i=%d, j=%d, F=%d, Fs=%d\n", i, j, F[i*m+j], Fs[i*m+j]);
				passed = 0;
			}
		}
	}
	if(passed)
		printf("result verified\n");
  return 0;
}
