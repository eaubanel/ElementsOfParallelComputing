/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm  3.5 from
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
 * Implementation of Conway's Game of Life, with periodic
 * boundary conditions
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void updateGrid(char **grid, char **newGrid, int n);
void display(char **grid, int n);
char **initialize(int seed, int n);
char **allocate2D(int n);

int main(int argc, char **argv){
	int n; // n x n grid
	int ngen; // number of generations
	char **grid; // Game grid
	char **newGrid; //copy of grid

	if(argc <2){
		fprintf(stderr,"usage: %s n [seed]\n", argv[0]);
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	int seed;
	if(3 == argc)
		seed = strtol(argv[2], NULL, 10);
	else
		seed = -1;
	grid = initialize(seed, n);
	newGrid = allocate2D(n);
	if(!grid || !newGrid){
		fprintf(stderr,"couldn't allocate memory for grid\n");
		return 1;
	}
	display(grid, n);
	printf("enter number of generations: ");
	scanf("%d", &ngen);
	getchar(); // consume newline

	while(ngen > 0){
		updateGrid(grid, newGrid, n);
		char ** temp = grid;
		grid = newGrid;
		newGrid = temp;
		display(grid, n);
		getchar(); // pause program until return pressed
		ngen--;
	}
	return 0;
}

void updateGrid(char **grid, char **newGrid, int n){
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++){
			int up = (i-1+n)%n;
			int down = (i+1)%n;
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

void display(char **grid, int n){
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			char alive = 'o';
			char dead = '.';
			printf("%c",grid[i][j]?alive:dead);
		}
		putchar('\n');
	}
}

char **initialize(int seed, int n){
	char *temp = calloc(n*n, sizeof(char));
	char **grid = malloc(n*sizeof(char *));
	if(!temp || !grid)
		return NULL;
	if(-1 == seed)
		srand(time(NULL));
	else
		srand(seed);
	for(int i=0; i<n*n; i++)
		if(rand() > RAND_MAX/2)
			temp[i] = 1;
	for(int i=0; i<n; i++)
		grid[i] = temp + i*n;
	return grid;
}
	
char **allocate2D(int n){
		char *temp = malloc(n*n*sizeof(char));
		char **a = malloc(n*sizeof(char *));
		if(!temp || !a)
			return NULL;
		for(int i=0; i<n; i++)
			a[i] = temp + i*n;
		return a;
}


