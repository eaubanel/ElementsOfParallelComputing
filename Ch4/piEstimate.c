/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.4 of
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
 * Implementation of recursive estimation of pi
 */
 #include <stdio.h>
#include <stdlib.h>
#include <time.h>

float recPi(int n);
int cutoff; //cutoff for recursion

int main(int argc, char **argv){
	struct timespec tstart,tend; 
  float timer;

	if(argc < 3){
		fprintf(stderr,"usage: %s n cutoff\n", argv[0]);
		return 1;
	}
	int n = strtol(argv[1], NULL, 10);
	cutoff = strtol(argv[2], NULL, 10);

	clock_gettime(CLOCK_MONOTONIC, &tstart);
	float piEst = recPi(n)*4/n;
	clock_gettime(CLOCK_MONOTONIC, &tend);
  timer = (tend.tv_sec-tstart.tv_sec) +
        (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;

	printf("pi is approx %f\n", piEst);
	printf("time in s: %f\n", timer);
	return 0;
}

float piEst(int n){
	float sum=0.0;
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	unsigned int seed = t.tv_nsec;
	printf("seed=%u\n", seed);
	for(int i=0; i<n; i++){
		float x = (float)rand_r(&seed)/RAND_MAX*2-1;
		float y = (float)rand_r(&seed)/RAND_MAX*2-1;
		if(x*x + y*y <= 1.0)
			sum++;
	}
	return sum;
}

float recPi(int n){
	float sum;
	if(n < cutoff)
		sum = piEst(n);
	else{
		float sum1 = recPi(n/2);
		float sum2 = recPi(n/2);
		sum = sum1 + sum2;
	}
	return sum;
}

