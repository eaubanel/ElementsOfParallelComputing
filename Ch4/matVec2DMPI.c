/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.18 from
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
 * MPI implementation of 2D decomposed nxn matrix-vector multiplication,
 * where p processors divide n, and p is square( = q^2)
 * Fixes bug in Algorithm 4.18, by properly defining process groups for
 * broadcast to redistribute vector, and by using two broadcasts.
 * Each process generates random block of matrix;
 * Vector generated randomly in such a way that it is scattered across
 * elements in a row, but replicated across all rows.
 * After multiplication, result vector is distributed in same fashion
 * as initial vector.
 * Result is verified by comparing result of sequential multiplication
 * with result of parallel multiplication gathered from first row of processes.
 * See Fig. 4.13-4.15 for diagrams for decomposition and communication.
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
// sequential multiplication for verification
// returns result vector
float *verifyMatVec(int n, int p);

int main(int argc, char **argv){
	float **A; //matrix
	float *b; //right-hand-side vector
	float *c; //local result vector
	float *cs; //result from sequential multiplication
	float *cp; // result gathered from parallel multiplication
	float *cr; //global result vector
	int id; //my id
	int p; //number of processes
	double time;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int q = sqrt(p);
	if(q*q != p){
		if(!id) fprintf(stderr,"p must be square\n");
		MPI_Finalize();
		return 1;
	}
	if(argc <2){
		if(!id) fprintf(stderr,"usage: %s n\n", argv[0]);
		MPI_Finalize();
		return 1;
	}
	int n = strtol(argv[1], NULL, 10);
	if(n%q){
		if(!id) fprintf(stderr,"n must be divisible by q\n");
		MPI_Finalize();
		return 1;
	}
	int nb = n/q; //dimension of block 
	//coordinates of each process in 2D grid
	int rowID = id/q;
	int colID = id%q;

	b = malloc(nb*sizeof(float));
	c = malloc(nb*sizeof(float));
	cr = malloc(nb*sizeof(float));
	A = allocate2D(nb, nb);
	if(!b || !c || !cr || !A){
		if(!id) fprintf(stderr,"couldn't allocate memory\n");
		MPI_Finalize();
		return 1;
	}

	srand(id*p);
	for(int i=0; i<nb; i++){
		for(int j=0; j<nb; j++){
			A[i][j] = rand();
		}
	}
	srand(colID*q);
	for(int i=0; i<nb; i++){
		b[i] = rand();
	}

	MPI_Barrier(MPI_COMM_WORLD);
	time = -MPI_Wtime();

	matvec(A, b, c, nb, nb);

	//reduce partial result c across each row to first column in array cr
	MPI_Comm rowComm;
	//communicator for each row; colID is rank in new communicator
	MPI_Comm_split(MPI_COMM_WORLD, rowID, colID, &rowComm);
	MPI_Reduce(c, cr, nb, MPI_FLOAT, MPI_SUM, 0, rowComm);
		
	//Redistribution of result vector
	MPI_Group worldGroup, newGroup;
	MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
	MPI_Comm tempComm, bcastComm;
	int ranks[q+1]; // processes in a group
	//First broadcast to columns 1 to q-1	
	for(int j=1; j<q; j++){
		ranks[0] = j*q;
		for(int i=0; i<q; i++)
			ranks[i+1] = j + i*q;
		MPI_Group_incl(worldGroup, q+1, ranks, &newGroup);
		MPI_Comm_create(MPI_COMM_WORLD, newGroup, &tempComm);
		if(j == colID || j*q == id){
			bcastComm = tempComm;
		}
	}
	if(id) //process 0 doesn't take part
		MPI_Bcast(cr, nb, MPI_FLOAT, 0, bcastComm);
	//Now broadcast in first column
	for(int i=0; i<q; i++)
			ranks[i] = i*q;
	MPI_Group_incl(worldGroup, q, ranks, &newGroup);
	MPI_Comm_create(MPI_COMM_WORLD, newGroup, &bcastComm);
	//only processes in first column take part
	if(0 == colID)
			MPI_Bcast(cr, nb, MPI_FLOAT, 0, bcastComm);

	time += MPI_Wtime();
	double ptime;
	MPI_Reduce(&time, &ptime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!id){
		printf("time in seconds: %f\n", ptime);
		//verification
 		cs = verifyMatVec(n, p);
		cp = malloc(n*sizeof(float));
	}
	if(!id && (!cp || !cs))
			fprintf(stderr, "can't allocate memory for sequential multiplication\n");
	else{
		// gather parallel result across first row; result in process 0
		if(0 == rowID)
			MPI_Gather(cr, nb, MPI_FLOAT, cp, nb, MPI_FLOAT, 0, rowComm);
		if(!id){
			printf("machine epsilon = %g\n", FLT_EPSILON);
			float max = 0.0;
			for(int i=0; i<n; i++){
				float err = fabs(cp[i]-cs[i])/cs[i];
				if(err > max)
					max = err;
			}
			printf("maximum relative difference: %g\n", max);
		}
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

float *verifyMatVec(int n, int p){
	float *c = malloc(n*sizeof(float));
	float *b = malloc(n*sizeof(float));
	float **A = allocate2D(n, n);
	if(!c || !b || !A){
		fprintf(stderr,"couldn't allocate memory for sequential verification\n");
		return NULL;
	}
	int q = sqrt(p);
	int nb = n/q;
	for(int id=0; id<p; id++){
		int rowID = id/q;
		int colID = id%q;
		srand(id*p);
		for(int i=0; i<nb; i++){
			for(int j=0; j<nb; j++){
				A[i+rowID*nb][j+colID*nb] = rand();
			}
		}
	}
	for(int colID=0; colID<q; colID++){
		srand(colID*q);
		for(int i=0; i<nb; i++)
			b[i+colID*nb] = rand();
	}
	
	matvec(A, b, c, n, n);
	free(A);
	free(b);
	return c;
}
