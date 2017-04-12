/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm  3.1 from
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
 * Implementation of 2-dimensional k-means algorithm
 * Reads a file of points (with x and y coords on each line)
 * specified in the command line and prompts for number of 
 * clusters and initial guess of each cluster center. 
 * Outputs assignment of each vector to a cluster.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct{
	double *x, *y;
} Vector;


// allocate struct of 2 arrays of doubles of length n
void allocVector(Vector *v, int n);
// Return the index of the closest cluster center to (x,y)
int findClosest(double x, double y, Vector cluster, int k);

int main(int argc, char **argv){
	int n; // number of vectors
	int k; // numbef of clusters 
	Vector vector; // vectors
	Vector cluster; // cluster centers
	Vector clusterNew; // for use in calculating new cluster centers
	int *clusterSize; // array of cluster sizes
	int *closest; //assignment of vectors to cluster centers

	if(argc <2){
		fprintf(stderr,"usage: %s filename\n", argv[0]);
		return 1;
	}
	FILE *f = fopen(argv[1], "r");
	if(f == NULL){
		fprintf(stderr,"can't open file %s\n", argv[1]);
		return 1;
	}
	//First line of file should indicate number of points
	fscanf(f, "%d", &n);

	//Read points and initialize closest
	allocVector(&vector, n);
	if((closest = malloc(n*sizeof(int))) == NULL){
		fprintf(stderr,"couldn't allocate array of %d ints\n", n);
		return 1;
	}
	for(int i=0; i<n; i++){
		if(fscanf(f, "%lf %lf", &vector.x[i], &vector.y[i]) != 2){
			fprintf(stderr,"something wrong with data in file\n");
			return 1;
		}
		closest[i] = 0;
	}
	
	// Read clusters and initialize cluster arrays
	printf("enter number of clusters: ");
	scanf("%d", &k);
	allocVector(&cluster, k);
	allocVector(&clusterNew, k);
	if((clusterSize = malloc(k*sizeof(int))) == NULL){
		fprintf(stderr,"couldn't allocate array of %d ints\n", k);
		return 1;
	}
	printf("enter coordinates for guess of %d clusters\n", k);
	for(int i=0; i<k; i++){
		if(scanf("%lf %lf", &cluster.x[i], &cluster.y[i]) != 2){
			fprintf(stderr,"something wrong with coordinate\n");
			return 1;
		}
		clusterNew.x[i] = clusterNew.y[i] = 0.0;
		clusterSize[i] = 0;
	}

	// k-means
	bool converged = false;
	while(!converged){
		converged = true;
		for(int j=0; j<n; j++){
			int i = findClosest(vector.x[j], vector.y[j], cluster, k);
			if(i != closest[j])
				converged = false;
			closest[j] = i;
			clusterNew.x[i] += vector.x[j];
			clusterNew.y[i] += vector.y[j];
			clusterSize[i]++;
		}
		for(int i=0; i<k; i++){
			cluster.x[i] = clusterNew.x[i]/clusterSize[i];
			cluster.y[i] = clusterNew.y[i]/clusterSize[i];
			clusterNew.x[i] = clusterNew.y[i] = 0.0;
			clusterSize[i] = 0.0;
		}
	}

	for(int j=0; j<n; j++)
		printf("%d\n", closest[j]);

	return 0;
}

void allocVector(Vector *v, int n){
	v->x = malloc(n*sizeof(double));
	v->y = malloc(n*sizeof(double));
	if(v->x == NULL || v->y == NULL){
		fprintf(stderr,"couldn't allocate memory for %d doubles\n", n);
		exit(1);
	}
}

int findClosest(double x, double y, Vector cluster, int k){
	double min = sqrt((x-cluster.x[0])*(x-cluster.x[0]) + (y-cluster.y[0])*(y-cluster.y[0]));
	int mini = 0;
	for(int i=1; i<k; i++){
		double d = sqrt((x-cluster.x[i])*(x-cluster.x[i]) + (y-cluster.y[i])*(y-cluster.y[i]));
		if(d < min){
			min = d;
			mini = i;
		}
	}
	return mini;
}
