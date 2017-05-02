/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.9 of
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
 * Remove duplicates of list of integers of limited range
 * using OpenMP and compare-and-swap (__sync_bool_compare_and_swap)
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>

int main(int argc, char **argv){
	int n; //number of values
	int R; //values in range [0..R)
	int *a; //array of values
	char *t; //array for marking duplicates
	char *ts; //array for verification of parallel loop
	struct timespec tstart,tend; 
  float timer;

	if(argc < 3){
		fprintf(stderr,"usage: %s n R\n", argv[0]);
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	R = strtol(argv[2], NULL, 10);
	a = malloc(n*sizeof(int));
	t = calloc(R, sizeof(char));
	ts = calloc(R, sizeof(char));
	if(a == NULL || t == NULL || ts == NULL){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	srand(time(NULL));
	for(int i=0; i<n; i++)
		a[i] = rand()%R;

	clock_gettime(CLOCK_MONOTONIC, &tstart);
	#pragma omp parallel for
	for(int i=0; i<n; i++)
		__sync_bool_compare_and_swap(t+a[i],0,1);
	clock_gettime(CLOCK_MONOTONIC, &tend);
  timer = (tend.tv_sec-tstart.tv_sec) +
        (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("time to mark duplicates in s: %f\n", timer);

	int k=0;
	for(int i=0; i<R; i++)
		if(t[i])
			a[k++] = i;
	printf("%d distinct values\n", k);

	//verification
	for(int i=0; i<n; i++)
		ts[a[i]] = 1;
	int passed = 1;
	for(int i=0; i<R; i++)
		if(t[i] != ts[i]){
			printf("i=%d, t[i]=%d, ts[i]=%d\n", i, t[i], ts[i]);
			passed = 0;
		}
	if(passed)
		printf("result verified\n");

	return 0;	
}

