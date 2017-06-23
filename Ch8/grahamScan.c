/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 8.3
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
 * Implementation of Graham Scan for 2D convex hull. Inputs point 
 * coordinates in format x,y for each point, where x and y are integers. 
 * Sorts points lexicographically, computes upper hull, then outputs
 * points on upper hull.
 */
#include <stdio.h>
#include <stdlib.h>

typedef struct{
	int num; //point number
	int x, y; //x,y coordinate of point
} PointT;

//area of triangle abc is negative if a-b-c forms a right hand turn
#define RIGHT(a, b, c) 	(((b).x - (a).x)*((c).y - (a).y) - \
												((c).x - (a).x)*((b).y-(a).y) < 0)

int comparePoints(const void *q1, const void *q2);

int main(int argc, char **argv){
	PointT *S; //array of points
	PointT *H; //array of points on hull
	int n; //number of points in set
	int m; //number of points on hull
	
	if(argc < 2){
		fprintf(stderr,"usage: %s n\n", argv[0]);
		return 1;
	}
	n = strtol(argv[1], NULL, 10);
	m = 0;
	S = malloc(n*sizeof(PointT));
	H = malloc(n*sizeof(PointT));
	if(!S || !H){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}

	for(int i=0; i<n; i++){
		if(scanf("%d,%d", &S[i].x, &S[i].y) != 2){
			fprintf(stderr, "error in reading %d points\n", n);
			return 1;
		}
		S[i].num = i;
	}

	qsort(S, n, sizeof(PointT), comparePoints);

	//S[0] is leftmost point, so is on hull
	H[0] = S[0];
	H[1] = S[1];
	m = 2;
	for(int i=2; i<n; i++){
		H[m++] = S[i];
		while(m > 2 && !RIGHT(H[m-3], H[m-2], H[m-1])){
			//remove second to last point in H
			H[m-2] = H[m-1];
			m--;
		}
	}

	for(int i=0; i<m; i++)
		printf("%d: %d,%d\n", H[i].num, H[i].x, H[i].y);

	return 0;
}

int comparePoints(const void *q1, const void *q2){
	PointT p1 = *((PointT *) q1);
	PointT p2 = *((PointT *) q2);
	if(p1.x != p2.x) 
		return p1.x - p2.x;
	else 
		return p1.y - p2.y;
}
