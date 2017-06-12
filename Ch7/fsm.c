/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 7.1
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
 * Sequential implementation of Fast Sweeping Method for solution of
 * two-dimensional Eikonal equation.
 * Reads speed function at each point in 2D grid, then reads 
 * coordinates of boundary of front, terminated with negative value.
 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MIN(x, y) ((x)<(y)?(x):(y))

//Sweep through grid in direction indicated by starting and ending indices
//given by ia, ib, ja, jb. Note the extra layer around grid to handle external
//boundary conditions. This affects indexing (e.g. left-to-right sweep in i 
//goes from 1 to ni). Also finds maximum relative difference in solution 
//(maxerr) after sweep.
void sweep(double **U, int ia, int ib, int ja, int jb, double **F, double h, double *maxerr);
//Initialize F and U 2D arrays
void initialize(double **F, double **U, int ni, int nj);
double **allocate2D(int nrows, int ncols);

int main(int argc, char **argv){
	int ni, nj; //number of rows and columns
	double h; //grid spacing
	double **F; //speed function at each point on grid
	double **U; //solution matrix
	const double tol = 1e-6; //tolerance for convergence
	double maxerr; //maximum relative error

	if(argc < 4){
		fprintf(stderr,"usage: %s nrows ncols grid_spacing\n", argv[0]);
		return 1;
	}
	ni = strtol(argv[1], NULL, 10);
	nj = strtol(argv[2], NULL, 10);
	h = strtod(argv[3], NULL);
	F = allocate2D(ni+2, nj+2);
	U = allocate2D(ni+2, nj+2);
	if(!F || !U){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	initialize(F, U, ni, nj);
	
	do{
		maxerr = 0.0;
		sweep(U, ni, 1, 1, nj, F, h, &maxerr);
		sweep(U, ni, 1, nj, 1, F, h, &maxerr);
		sweep(U, 1, ni, nj, 1, F, h, &maxerr);
		sweep(U, 1, ni, 1, nj, F, h, &maxerr);
	} while(maxerr > tol);
	
	for(int i=1; i<=ni; i++){
		for(int j=1; j<=nj; j++){
			printf("%g ", U[i][j]);
		}
		putchar('\n');
	}

	return 0;
}

inline static double solveQuadratic(double **U, int i, int j, double **F, double h){
	if(0 == U[i][j]) return 0;
	double a = MIN(U[i-1][j], U[i+1][j]);
	double b = MIN(U[i][j-1], U[i][j+1]);
	double hinvF = h/F[i][j];
	double amb = a - b;
	if(fabs(amb) >= hinvF)
		return MIN(a, b) + hinvF;
	else
		return (a + b + sqrt(2*(hinvF*hinvF) - amb*amb))/2.0;
}

void sweep(double **U, int ia, int ib, int ja, int jb, double **F, double h, double *maxerr){
	int stepi = ib<ia?-1:1;
	int stepj = jb<ja?-1:1;
	for(int i=ia; stepi*i<=stepi*ib; i+=stepi){
		for(int j=ja; stepj*j<=stepj*jb; j+=stepj){
			double unew = solveQuadratic(U, i, j, F, h);
			if(unew < U[i][j]){
				double err = fabs(U[i][j]-unew)/U[i][j];
				if(err > *maxerr) *maxerr = err;
				U[i][j] = unew;
			}
		}
	}
}

void initialize(double **F, double **U, int ni, int nj){
	for(int i=0; i<ni+2; i++){
		for(int j=0; j<nj+2; j++){
			if(scanf("%lf", &F[i][j]) != 1){
				fprintf(stderr, "invalid input (F)\n");
				exit(1);
			}
			U[i][j] = DBL_MAX;
		}
	}
	int i, j;
	while(1){
		int sr = scanf("%d %d", &i, &j);
		//i = -1 is sentinel
		if(i < 0) break;
		if(sr != 2){
				fprintf(stderr, "invalid input (U)\n");
				exit(1);
		}
		U[i][j] = 0.0;
	}
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
