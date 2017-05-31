/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm  5.3 from
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
 * Parallel MPI implementation of Conway's Game of Life, with periodic
 * boundary conditions, and extra layer of ghost cells.
 * Also includes initialization and display of grid.
 * Assumes number of rows (n) divisible by number of processors (p).
 * Process 0 has entire grid, plus two ghost rows (first row).
 * Other processes have grid with n/p + 4 rows (includes 4 ghost rows).
 * If DISPLAY defined, displays grid after every DISP_FREQ generations,
 * otherwise only displays initial and final grids.
 * Grid diplayed by gathering it from all processes, then 
 * process 0 writes it to file gameOfLifeMPI.txt.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

//update grid, store results in newGrid, for m+2 rows (start=0) or 
//m rows (start=1), and n columns
void updateGrid(char **grid, char **newGrid, int start, int m, int n);
//writes entire grid to file f
void display(FILE *f, char **grid, int n);
//called by process 0 to initialize grid with random values
//using seed for rand(), or time(NULL) if seed is -1
char **initialize(int seed, int n);
//allocate char array with m rows and n columns
char **allocate2D(int m, int n);

#define DISP_FREQ 10

int main(int argc, char **argv){
	int n; // n x n grid
	int ngen; // number of generations
	char **grid; // Game grid
	char **newGrid; //copy of grid
	int id; //my id
	int p; //number of processes
	FILE *f; //grids written to this file
	double time;

	MPI_Status status;
  MPI_Request req_send_up, req_send_down;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if(argc <2){
		if(!id) fprintf(stderr,"usage: %s n [seed]\n", argv[0]);
		MPI_Finalize();
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	if(n%p){
		if(!id) fprintf(stderr,"n must be divisible by number of processes\n");
		MPI_Finalize();
		return 1;
	}
	int m = n/p; //number of rows per process
	
	int seed;
	if(!id){
		if(3 == argc)
			seed = strtol(argv[2], NULL, 10);
		else
			seed = -1;
		grid = initialize(seed, n);
	} else{
			grid = allocate2D(m+4, n);
	}
	if(grid == NULL){
			if(!id) fprintf(stderr,"couldn't allocate memory for grid\n");
			MPI_Finalize();
	}

	MPI_Scatter(&grid[2][0], m*n, MPI_CHAR, &grid[2][0], m*n, MPI_CHAR,
							0, MPI_COMM_WORLD);

	if(!id)
		newGrid = allocate2D(n+2, n);
	else
		newGrid = allocate2D(m+4, n);
	if(newGrid == NULL){
		if(!id) fprintf(stderr,"couldn't allocate memory for newGrid\n");
		MPI_Finalize();
		return 1;
	}
	if(!id){
		f = fopen("gameOfLifeMPI.txt","w");
		display(f, grid, n);
		printf("enter number of generations: ");
		fflush(stdout);
		scanf("%d", &ngen);
	}
	MPI_Bcast(&ngen, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int nbDown = (id+1)%p; //lower neighbour
	int nbUp = (id-1+p)%p; //upper neighbour
	MPI_Barrier(MPI_COMM_WORLD);
	time = -MPI_Wtime();
	for(int k=0; k<ngen; k++){
		int offset = k%2;
		if(!offset){
			//exhange boundary values
			MPI_Isend(&grid[m][0], 2*n, MPI_CHAR, nbDown, 1, 
							MPI_COMM_WORLD, &req_send_down);
			MPI_Isend(&grid[2][0], 2*n, MPI_CHAR, nbUp, 2, 
							MPI_COMM_WORLD, &req_send_up);
			MPI_Recv(&grid[m+2][0], 2*n, MPI_CHAR, nbDown, 2, 
							MPI_COMM_WORLD, &status);
			MPI_Recv(&grid[0][0], 2*n, MPI_CHAR, nbUp, 1, 
							MPI_COMM_WORLD, &status);
		}
		updateGrid(grid, newGrid, offset, m, n);
		char ** temp = grid;
		grid = newGrid;
		newGrid = temp;

		#ifdef DISPLAY
		if(!id && k % DISP_FREQ == 0){
			MPI_Gather(&grid[2][0], m*n, MPI_CHAR, &grid[2][0], m*n, MPI_CHAR,
								0, MPI_COMM_WORLD);
		
			display(f, grid, n);
		}
		#endif
	}
	time += MPI_Wtime();
	double ptime;
	MPI_Reduce(&time, &ptime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!id)
		printf("time in seconds: %f\n\n", ptime);
	MPI_Gather(&grid[2][0], m*n, MPI_CHAR, &grid[2][0], m*n, MPI_CHAR,
						0, MPI_COMM_WORLD);
	if(!id)
		display(f, grid, n);
	MPI_Finalize();
	return 0;
}

void updateGrid(char **grid, char **newGrid, int start, int m, int n){
	for(int i=1+start; i<=m+2-start; i++)
		for(int j=0; j<n; j++){
			int up = i-1;
			int down = i+1;
			int left = (j-1+n)%n;
			int right = (j+1)%n;
			int sumAlive = grid[up][left] + grid[up][j] + grid[up][right] +
					grid[i][left] + grid[i][right] + 
					grid[down][left] + grid[down][j] + grid[down][right];
			if(0 == grid[i][j] && 3 == sumAlive)
				newGrid[i][j] = 1;
			else if( 1 == grid[i][j] && (2 == sumAlive || 3 == sumAlive))
				newGrid[i][j] = 1;
			else
				newGrid[i][j] = 0;
		}
}

void display(FILE *f, char **grid, int n){
	for(int i=2; i<n+2; i++){
		for(int j=0; j<n; j++){
			char alive = 'o';
			char dead = '.';
			fprintf(f, "%c",grid[i][j]?alive:dead);
		}
		fputc('\n', f);
	}
}

char **initialize(int seed, int n){
	int nrows = 2+n;
	char *temp = calloc(nrows*n, sizeof(char));
	char **grid = malloc(nrows*sizeof(char *));
	if(!temp || !grid)
		return NULL;
	if(-1 == seed)
		srand(time(NULL));
	else
		srand(seed);
	//leave first two rows blank (ghost rows)
	for(int i=2*n; i<nrows*n; i++)
		if(rand() > RAND_MAX/2)
			temp[i] = 1;
	for(int i=0; i<nrows; i++)
		grid[i] = temp + i*n;
	return grid;
}
	
char **allocate2D(int m, int n){
		char *temp = malloc(m*n*sizeof(char));
		char **a = malloc(m*sizeof(char *));
		if(!temp || !a)
			return NULL;
		for(int i=0; i<m; i++)
			a[i] = temp + i*n;
		return a;
}


