/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.17 from
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
 * MPI implementation of row-wise nxn matrix-vector multiplication
 * Process 0 generates random matrix then scatters to all processors
 * Process 0 broadcasts vector to all processes.
 * After multiplication, result vector is gathered at process 0
 * then broadcast to all processes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "mpi.h"

//allocate char array with m rows and n columns
float **allocate2D(int m, int n);
// Multiply m rows of matrix A by vector b of length n
// result in vector c
void matvec(float **A, float *b, float *c, int m, int n);

int main(int argc, char **argv){
	//process 0 has nxn matrix, others have n*n/p matrix
	float **A; //matrix
	float *b; //right-hand-side vector
	float *c; //result vector
	float *cs; //result vector from sequential multiplication
	int id; //my id
	int p; //number of processes
	double time;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if(argc <2){
		if(!id) fprintf(stderr,"usage: %s n\n", argv[0]);
		MPI_Finalize();
		return 1;
	}
	int n = strtol(argv[1], NULL, 10);
	if(n%p){
		if(!id) fprintf(stderr,"n must be divisible by number of processes\n");
		MPI_Finalize();
		return 1;
	}
	int m = n/p; //number of rows per process

	b = malloc(n*sizeof(float));
	c = malloc(n*sizeof(float));
	if(!id){
		cs = malloc(n*sizeof(float));
		A = allocate2D(n, n);
	} else{
		A = allocate2D(m, n);
	}
	if(!b || !c || !A){
		if(!id) fprintf(stderr,"couldn't allocate memory\n");
		MPI_Finalize();
		return 1;
	}

	if(!id){
		for(int i=0; i<n; i++){
			b[i] = rand();
				for(int j=0; j<n; j++){
					A[i][j] = rand();
				}
		}
		//sequential multiplication (for verification)
		matvec(A, b, cs, n, n);
	}

	MPI_Scatter(&A[0][0], m*n, MPI_FLOAT, &A[0][0], m*n, MPI_FLOAT,
						0, MPI_COMM_WORLD);
	MPI_Bcast(b, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	time = -MPI_Wtime();

	matvec(A, b, c, m, n);

	//these two calls could be replaced by all-gather
	MPI_Gather(c, m, MPI_FLOAT, c, m, MPI_FLOAT, 
					0, MPI_COMM_WORLD);
	MPI_Bcast(c, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	time += MPI_Wtime();
	double ptime;
	MPI_Reduce(&time, &ptime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!id){
		printf("time in seconds: %f\n\n", ptime);
		printf("machine epsilon = %g\n", FLT_EPSILON);
		float max = 0.0;
		for(int j=0;j<n;j++){
			float err = fabs(c[j]-cs[j])/cs[j];
			if(err > max)
				max = err;
		}
		printf("maximum relative difference: %g\n", max);
	}
	MPI_Finalize();
  return 0;
}

void matvec(float **A, float *b, float *c, int m, int n){
	for(int i=0; i<m; i++){
		c[i] = 0;
		for(int j=0; j<n; j++){
			c[i] += A[i][j]* b[j];
		}
	}
}

float **allocate2D(int m, int n){
		float *temp = malloc(m*n*sizeof(float));
		float **a = malloc(m*sizeof(float *));
		if(!temp || !a)
			return NULL;
		for(int i=0; i<m; i++)
			a[i] = temp + i*n;
		return a;
}
