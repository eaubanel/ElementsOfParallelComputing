/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 8.2
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
 * Implementation of Merge Hull for 2D convex hull. Inputs point 
 * coordinates in format x,y for each point, where x and y are integers. 
 * Sorts points lexicographically, computes hull, then outputs
 * points on hull.
 */
#include <stdio.h>
#include <stdlib.h>

typedef struct{
	int num; //point number
	int x, y; //x,y coordinate of point
} PointT;

typedef struct{
	PointT *H; //pointer to hull array
	int nu, nl; //number of points on upper and lower hull
} HullT;

#define AREA(a, b, c) 	(((b).x - (a).x)*((c).y - (a).y) - \
												((c).x - (a).x)*((b).y-(a).y))
//area of triangle abc is negative if a-b-c forms a right hand turn
#define RIGHT(a, b, c) 	(AREA(a, b, c) < 0)
#define LEFT(a, b, c) 	(AREA(a, b, c) > 0)

//Computes hull of n points S, returns hull.
//Alternatively merges between arrays H and Hbuff.
//Hull ordered clockwise starting from leftmost point.
HullT mergeHull(PointT *S, int n, PointT *H, PointT *Hbuff);
//Join two hulls HI1 (left) and HI2 (right), storing new hull 
//in HBuff and returning hull.
HullT joinHulls(HullT HI1, HullT HI2, PointT *Hbuff);
//Find tangent between point p on left and right 
// upper hull (upper=1) or lower hull (upper=0) of HI
PointT *findTangent(PointT *p, HullT HI, int upper);
//Needed for qsort
int comparePoints(const void *q1, const void *q2);

int main(int argc, char **argv){
	PointT *S; //array of points
	PointT *H; //array of points on hull
	PointT *Hbuff; //copy of array of points on hull
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
	Hbuff = malloc(n*sizeof(PointT));
	if(!S || !H || !Hbuff){
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

	HullT HI = mergeHull(S, n, H, Hbuff);
	H = HI.H;

	for(int i=0; i<HI.nu+HI.nl; i++)
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

HullT mergeHull(PointT *S, int n, PointT *H, PointT *Hbuff){
	//base case has 2 or 3 points, which are stored in clockwise order
	if(n <= 3){
		HullT HI;
		HI.H = H;
		HI.H[0] = S[0];
		if(2 == n){
			HI.nu = HI.nl = 1;
			HI.H[1] = S[1];
		}else{
			if(RIGHT(S[0], S[1], S[2])){
				HI.nu = 2;
				HI.nl = 1;
				HI.H[1] = S[1];
				HI.H[2] = S[2];
			}else if(LEFT(S[0], S[1], S[2])){
				HI.nu = 1;
				HI.nl = 2;
				HI.H[1] = S[2];
				HI.H[2] = S[1];
			}else{ //three points are colinear; remove middle point
				HI.nu = HI.nl = 1;
				HI.H[1] = S[2];
			}
		}
		return HI;
	}else{
		int mid = n/2;
		HullT HI1 = mergeHull(S, mid, Hbuff, H);
		int n1 = HI1.nu+HI1.nl;
		HullT HI2 = mergeHull(S+mid, n-mid, Hbuff+n1, H+n1);
		return joinHulls(HI1, HI2, H);
	}
}

HullT joinHulls(HullT HI1, HullT HI2, PointT *Hbuff){
	PointT *q;
	PointT *H1 = HI1.H;
	int n1 = HI1.nu + HI1.nl; //number of points on left hull
	//find upper tangent (pu, qu)
	int min = 0;
	int max = HI1.nu; //search ends at rightmost point
	int mid;
	while(min <= max){
		mid = (min+max)/2;
		q = findTangent(H1+mid, HI2, 1);//find upper tangent of H1[mid] and H2
		//check left
		if(mid > 0 && !RIGHT(H1[mid-1], H1[mid], *q))
			max = mid - 1;
		else if(mid == HI1.nu) //tangent found with rightmost point of H1
			break; 
		//check right
		else if(RIGHT(H1[mid], H1[mid+1], *q))
			min = mid + 1;
		else 
			break; //tangent found
	}
	PointT *pu = H1+mid;
	PointT *qu = q;
	//find lower tangent (pl, ql)
	min = HI1.nu;
	max = n1; //wrap around to first point
	while(min <= max){
		mid = (min+max)/2;
		q = findTangent(H1+mid%n1, HI2, 0);//find lower tangent of H1[mid] and H2
		//check right
		if(mid > HI1.nu && LEFT(H1[mid%n1], H1[mid-1], *q)) 
				max = mid - 1;
		else if(mid == n1) //tangent found with leftmost point of H1
			break;
		//check left
		else if(!LEFT(H1[(mid+1)%n1], H1[mid], *q))
			min = mid + 1;
		else 
			break; //tangent found
	}
	PointT *pl = H1+mid%n1;
	PointT *ql = q;
	//join hulls
	int m = 0;
	HullT HI;
	//upper hull from H1
	for(PointT *t=H1; t<=pu; t++)
		Hbuff[m++] = *t;
	PointT *qr = HI2.H + HI2.nu; //rightmost point on H2
	//upper hull from H2
	for(PointT *t=qu; t<qr; t++)
		Hbuff[m++] = *t;	
	HI.nu = m;
	//lower hull from H2
	if(ql == HI2.H){ //lower tangent on leftmost point of H2
		for(PointT *t=qr; t<HI2.H+HI2.nu+HI2.nl; t++)
			Hbuff[m++] = *t;
		Hbuff[m++] = *HI2.H;
	}else
		for(PointT *t=qr; t<=ql; t++)
			Hbuff[m++] = *t;
	//lower hull from H1
	if(pl != H1) //lower tangent not on leftmost point of H1
		for(PointT *t=pl; t<H1+n1; t++)
			Hbuff[m++] = *t;
	HI.nl = m - HI.nu;
	HI.H = Hbuff;
	return HI;
}

PointT *findTangent(PointT *p, HullT HI, int upper){
	PointT *H = HI.H;
	int n = HI.nu + HI.nl;
	int mid;
	if(upper){
		int min = 0;
		int max = HI.nu; //search ends at rightmost point
		while(min <= max){
			mid = (min+max)/2;
			//check left
			if(mid > 0 && RIGHT(*p, H[mid-1], H[mid]))
				max = mid - 1;
			else if(mid == HI.nu) //tangent found with rightmost point of hull
				break; 
			//check right
			else if(!RIGHT(*p, H[mid], H[mid+1]))
				min = mid + 1;
		 	else
				break; //tangent found	
		}
	}else{ //lower hull
		int min = HI.nu;
		int max = n; //wrap around to first point
		while(min <= max){
			mid = (min+max)/2;
			//check right
			if(mid > HI.nu && !LEFT(*p, H[mid%n], H[mid-1])) 
				max = mid - 1;
			else if(mid == n) //tangent found with leftmost point of hull
				break;
			//check left
			else if(LEFT(*p, H[(mid+1)%n], H[mid]))
				min = mid + 1;
			else
				break; //tangent found
		}
	}
	return H + mid%n;
}
