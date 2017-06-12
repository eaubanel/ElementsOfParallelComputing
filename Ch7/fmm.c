/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 7.2-7.3
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
 * Sequential implementation of Fast Marching Method for solution of
 * two-dimensional Eikonal equation.
 * Requires index min priority queue module (indexedMinPQ.h), which can
 * be adapted from indexed max priority queue (as a binary heap)
 * in Sedgewick's Algorithms in C, 3rd edition.
 * Reads speed function at each point in 2D grid, then reads 
 * coordinates of boundary of front, terminated with negative value.
 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "indexedMinPQ.h"

#define MIN(x, y) ((x)<(y)?(x):(y))
#define NOT_IN_DOMAIN(i, j, ni, nj) ((i)<1 || (i)>(ni) || (j)<1 || (j)>(nj))
//FAR not needed here but useful for parallel FMM
typedef enum {KNOWN, BAND, FAR} PointLabel;

//Initialize F, U, G 2D arrays and B array.
//Returns number of boundary points
int initialize(double **F, double **U, PointLabel **G, int *B, int nrows, int ncols);
//Update neighbours of point (i,j).
//Note the extra layer around grid to handle external
//boundary conditions. This affects indexing (e.g. visited elements in i direction
//are indexed from 1 to ni).
void updateNeighbors(double **U, PointLabel **G, int i, int j, int ni, int nj, double **F, double h);
double **allocate2D(int nrows, int ncols);
PointLabel **allocate2DLabel(int nrows, int ncols);

int main(int argc, char **argv){
	int ni, nj; //number of rows and columns
	double h; //grid spacing
	double L; //desired width of solution
	double **F; //speed function at each point on grid
	double **U; //solution matrix
	int *B; //record position of boundary points
	PointLabel **G; //label matrix (KNOWN, BAND, FAR)

	if(argc < 4){
		fprintf(stderr,"usage: %s nrows ncols grid_spacing width\n", argv[0]);
		return 1;
	}
	ni = strtol(argv[1], NULL, 10);
	nj = strtol(argv[2], NULL, 10);
	h = strtod(argv[3], NULL);
	L = strtod(argv[4], NULL);
	F = allocate2D(ni+2, nj+2);
	U = allocate2D(ni+2, nj+2);
	G = allocate2DLabel(ni+2, nj+2);
	B = malloc(ni*nj*sizeof(int));
	if(!F || !U || !G || !B){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	int nrows = ni+2, ncols = nj+2;
	int nB = initialize(F, U, G, B, nrows, ncols);
	//assume heap has no more than 2*(ni+nj) elements
	minHeapInit(U[0], nrows*ncols, 2*(ni+nj));

	for(int k=0; k<nB; k++)
		updateNeighbors(U, G, B[k]/ncols, B[k]%ncols, ni, nj, F, h);

	while(!minHeapEmpty()){
		int k = heapExtractMin();
		int i = k/ncols;
		int j = k%ncols;
		if(U[i][j] > L)
			break;
		G[i][j] = KNOWN;
		updateNeighbors(U, G, i, j, ni, nj, F, h); 
	}

	for(int i=1; i<=ni; i++){
		for(int j=1; j<=nj; j++){
			printf("%g ", U[i][j]);
		}
		putchar('\n');
	}

	return 0;
}

inline static double selectMin(PointLabel **G, double **U, int l, int m, int p, int q){
	double x = -1.0;
	if(KNOWN == G[l][m] && KNOWN == G[p][q])
		x = MIN(U[l][m], U[p][q]);
	else if(KNOWN != G[l][m] && KNOWN == G[p][q])
		x = U[p][q];
	else if(KNOWN == G[l][m] && KNOWN != G[p][q])
		x = U[l][m];
	return x;
}

inline static double solveQuadratic(PointLabel **G, double **U, int i, int j, double **F, double h){
	double a = selectMin(G, U, i+1, j, i-1, j);
	double b = selectMin(G, U, i, j+1, i, j-1);
	double hinvF = h/F[i][j];
	double amb = a-b;
	if(-1.0 == a)
		return b + hinvF;
	else if(-1.0 == b)
		return a + hinvF;	
	else if(fabs(amb) >= hinvF)
		return MIN(a, b) + hinvF;
	else
		return (a + b + sqrt(2*(hinvF*hinvF) - amb*amb))/2.0;
}

void updateNeighbors(double **U, PointLabel **G, int i, int j, int ni, int nj, double **F, double h){
	for(int l=i-1; l<=i+1; l++){
		for(int m=j-1; m<=j+1; m++){
			//only consider (i-1,j), (i+1,j), (i,j-1), (i,j+1)
			if((l != i && m != j) || (l == i && m == j)) continue;
			if((G[l][m] == KNOWN) || NOT_IN_DOMAIN(l,m,ni,nj)) continue;
			double utemp = solveQuadratic(G, U, l, m, F, h);
			if(utemp < U[l][m]){
				U[l][m] = utemp;
				G[l][m] = BAND;
				int k = l*(nj+2)+m;
				if(isInHeap(k))
					minHeapChange(k);
				else
					minHeapInsert(k);
			}
		}
	}
}

int initialize(double **F, double **U, PointLabel **G, int *B, int nrows, int ncols){
	for(int i=0; i<nrows; i++){
		for(int j=0; j<ncols; j++){
			if(scanf("%lf", &F[i][j]) != 1){
				fprintf(stderr, "invalid input (F)\n");
				exit(1);
			}
			U[i][j] = DBL_MAX;
			G[i][j] = FAR;
		}
	}
	int i=0, j=0, k=0;
	while(1){
		int sr = scanf("%d %d", &i, &j);
		//i = -1 is sentinel
		if(i < 0) break;
		if(sr != 2){
				fprintf(stderr, "invalid input (U)\n");
				exit(1);
		}
		U[i][j] = 0.0;
		G[i][j] = KNOWN;
		B[k++] = i*ncols + j;
	}
	return k;
}

double **allocate2D(int nrows, int ncols){
		double *temp = malloc(nrows*ncols*sizeof(double));
		double **a = malloc(nrows*sizeof(double *));
		if(!temp || !a)
			return NULL;
		for(int i=0; i<nrows; i++)
			a[i] = temp + i*ncols;
		return a;
}

PointLabel **allocate2DLabel(int nrows, int ncols){
		PointLabel *temp = malloc(nrows*ncols*sizeof(PointLabel));
		PointLabel **a = malloc(nrows*sizeof(PointLabel *));
		if(!temp || !a)
			return NULL;
		for(int i=0; i<nrows; i++)
			a[i] = temp + i*ncols;
		return a;
}
