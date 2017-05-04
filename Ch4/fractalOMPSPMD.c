/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.12 from
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
 * Generalized fractal, parallelized with OpenMP, SPMD round-robin style
 * Based on Gujar and Bhavsar, Computers & Graphics, 15(3):441-449, 1991.
 * Times execution if TIME defined (compile with -DTIME), 
 * otherwise not timed, and output written in PGM format
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <omp.h>

int main(int argc, char **argv){
  int niter=255; //maximum allowable number of iterations
  int threshold=10.0; //limit of |z|, beyond which z is said to diverge
  float len=3.0; //len^2 is area of picture
  float xmin=-1.5, ymin=-1.5;
  float ymax = ymin+len;
	char *count; //image stored in 1D array
	struct timespec tstart,tend; 
  float time;
  
	if(argc < 4){
		fprintf(stderr,"usage: %s n alpha chunk\n", argv[0]);
		return 1;
	}
  //number of points per column (and row) of image
	int n = strtol(argv[1], NULL, 10);
	//z = z^alpha + c
	float alpha = strtof(argv[2], NULL);
	//chunk size for round-robin scheduling
	int chunk = strtol(argv[3], NULL, 10);
  float ax = len/n;

  count = malloc(n*n*sizeof(char));
	if(count == NULL){
		fprintf(stderr,"couldn't allocate array of %d chars\n", n);
		return 1;
	}
#ifdef TIME
  clock_gettime(CLOCK_MONOTONIC, &tstart);
#endif
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
  	int nt = omp_get_num_threads();
  	int istart = id*chunk;
  	int iend = (id+1)*chunk-1;
		while(istart<n){
  		for(int i=istart;i<=iend;i++){
    		float cx = ax*i+xmin;
    		for(int j=0;j<n;j++){
      		float cy=ymax-ax*j;
      		float complex c=cx+I*cy;
					float complex z;
      		if(alpha>0)
        		z=0.0;
      		else
						z=1.0+I;
					int k;
      		for(k=1;k<=niter;k++)
        		if(cabsf(z)<threshold)
	  					z=cpowf(z,alpha)+c;
						else {
          		count[i*n+j] = k-1;
          		break;
      			}
      		if(k==niter+1) 
						count[i*n+j] = niter;
    		}
  		}
			istart += nt*chunk;
    	iend = istart + chunk - 1;
    	if(iend>n-1) 
				iend = n-1;
		}
	}
#ifdef TIME
  clock_gettime(CLOCK_MONOTONIC, &tend);
  time = (tend.tv_sec-tstart.tv_sec) +
        (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
  printf("time in s: %f\n", time);
#endif
#ifndef TIME
  printf("P2\n");
  printf("%d %d\n", n,n);
  printf("%d\n",niter);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++)
      printf("%d ",count[i*n+j]);
    printf("\n");
  }
#endif
}
