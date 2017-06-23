/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 8.1
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
 * Implementation of QuickHUll for 2D convex hull. Inputs point 
 * coordinates in format x,y for each point, where x and y are integers. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
	int num; //point number
	int x, y; //x,y coordinate of point
} PointT;

#define MIN(x, y) ((x)<(y)?(x):(y))
#define AREA(a, b, c) 	(((b).x - (a).x)*((c).y - (a).y) - \
												((c).x - (a).x)*((b).y-(a).y))
//area of triangle abc is negative if a-b-c forms a right hand turn
#define RIGHT(a, b, c) 	(AREA(a, b, c) < 0)
#define LEFT(a, b, c) 	(AREA(a, b, c) > 0)
//distance between point r and line p,q
#define DIST(r, p, q) (abs((q.x-p.x)*(p.y-r.y) - ((p.x-r.x)*(q.y-p.y))) \
											/sqrt((q.x-p.x)*(q.x-p.x) + (q.y-p.y)*(q.y-p.y)))

//Compute subhull for points above line p,q. Subhull also includes point p.
//Subhull stored in H. Size of subhull is returned in m.
//n points in array S partitioned in two and copied to array Sbuff
//During recursion partioning happens alternatively between S and Sbuff
//Array M used to mark points during partitioning
void subHull(PointT *S, PointT *Sbuff, int n, PointT p, PointT q, int *m, char *M, PointT *H);
//Partition n points in array S into two sets (both in Sbuff; pointer to second set returned),
//starting with points above p,q then those on the line or below.
//Size of each set returned in n1 and n2. Array M used to mark points during partitioning.
PointT *partition(PointT *S, PointT *Sbuff, int n, PointT p, PointT q, int *n1, int *n2, char *M);
//Return point with minimum x coordinate (breaks ties using min y coordinate).
PointT minX(PointT *S, int n);
//Return point with maximum x coordinate (breaks ties using max y coordinate).
PointT maxX(PointT *S, int n);
//Return point furthest away from line p,q.
PointT maxD(PointT *S, int n, PointT p, PointT q);

int main(int argc, char **argv){
	PointT *S; //array of points
	PointT *Sbuff; //working array to store points
	char *M; //array for marking points
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
	Sbuff = malloc(n*sizeof(PointT));
	M = malloc(n*sizeof(char));
	H = malloc(n*sizeof(PointT));
	if(!S || !Sbuff || !M || !H){
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

	PointT p = minX(S, n);
	PointT q = maxX(S, n);	
	int n1, n2;
	PointT *S1 = Sbuff;
	PointT *S2 = partition(S, Sbuff, n, p, q, &n1, &n2, M);
	int m1, m2; //size of upper and lower hulls
	//upper hull
	subHull(S1, Sbuff, n1, p, q, &m1, M, H);
	//lower hull
	subHull(S2, Sbuff, n2, q, p, &m2, M, H+m1);

	for(int i=0; i<m1+m2; i++)
		printf("%d: %d,%d\n", H[i].num, H[i].x, H[i].y);

	return 0;
}

void subHull(PointT *S, PointT *Sbuff, int n, PointT p, PointT q, int *m, char *M, PointT *H){
	if(0 == n){
		*m = 1;
		H[0] = p;
		return;
	}
	PointT r = maxD(S, n, p, q);
	int n1 = 0, n2 = 0;
	for(int i=0; i<n; i++){
		if(LEFT(p, r, S[i])){
			M[i] = 1;
			n1++;
		}else if(LEFT(r, q, S[i])){
			M[i] = 2;
			n2++;
		}else
			M[i] = 0;
	}
	PointT *S1 = Sbuff;
	PointT *S2 = Sbuff + n1;
	int i1 = 0, i2 = 0;
	for(int i=0; i<n; i++){
		if(1 == M[i])
			S1[i1++] = S[i];
		else if(2 == M[i])
			S2[i2++] = S[i];
	}
	int m1, m2;
	subHull(S1, S, n1, p, r, &m1, M, H); 
	subHull(S2, S, n2, r, q, &m2, M, H+m1); 
	*m = m1 + m2;
	return;
}

PointT *partition(PointT *S, PointT *Sbuff, int n, PointT p, PointT q, int *n1, int *n2, char *M){
	int nu = 0, nl = 0;
	for(int i=0; i<n; i++){
		if(LEFT(p, q, S[i])){
			nu++;
			M[i] = 1;
		}else{
			nl++;
			M[i] = 0;
		}
	}
	PointT *Su = Sbuff;
	PointT *Sl = Sbuff + nu;
	int iu = 0, il = 0;
	for(int i=0; i<n; i++){
		if(M[i])
			Su[iu++] = S[i];
		else
			Sl[il++] = S[i];
	}
	*n1 = nu;
	*n2 = nl;
	return Sl;
}

PointT minX(PointT *S, int n){
	PointT min = S[0];
	for(int i=1; i<n; i++)
		if(S[i].x < min.x)
		 	min = S[i];
		else if(S[i].x == min.x && S[i].y < min.y)
			min = S[i];	
	return min;
}

PointT maxX(PointT *S, int n){
	PointT max = S[0];
	for(int i=1; i<n; i++)
		if(S[i].x > max.x)
		 	max = S[i];
		else if(S[i].x == max.x && S[i].y > max.y)
			max = S[i];	
	return max;
}

PointT maxD(PointT *S, int n, PointT p, PointT q){
	PointT maxp = S[0];
	double max = DIST(S[0], p, q);
	for(int i=1; i<n; i++){
		double d = DIST(S[i], p, q);
		if(d > max){
			max = d;
			maxp = S[i];
		}
	}
	return maxp;
}

