/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 6.1
 * using indexed min priority queue, i.e. Dijkstra's algorithm, in
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
 * Sequential implementation of SSSP solver using Dijkstra's algorithm. 
 * Assumes nonzero integer weights.
 * Requires index min priority queue module (indexedMinPQ.h), which can
 * be adapted from indexed max priority queue (as a binary heap)
 * in Sedgewick's Algorithms in C, 3rd edition.
 * Reads graph as a list of weighted edges preceeded by 
 * two integers indicating the number of vertices and the number of edges.
 * Outputs list of distances to each vertex.
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "indexedMinPQ.h"

//read edge list from stdin and store CSR representation of graph
//in arrays V, E, and W.
void readGraph(int *V, int *E, int *W, int n, int m);

int main(int argc, char **argv){
	int n, m; //number of vertices and edges
	int *V, *E, *W; //Offset , edge, and weight arrays for CSR representation
	int s; //source vertex
	int *D; //distance array

	if(argc < 2){
		fprintf(stderr,"usage: %s source_vertex\n", argv[0]);
		return 1;
	}
	s = strtol(argv[1], NULL, 10);
  scanf("%d %d",&n, &m);
	if(s >= n || s < 0){
	  printf("invalid source vertex\n");
	  return 1;
  }
	V = malloc((n+1)*sizeof(int));
	E = malloc(m*sizeof(int));
	W = malloc(m*sizeof(int));
	D = malloc(n*sizeof(int));
	if(!V || !E || !W || !D){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	readGraph(V, E, W, n, m);

	for(int i=0; i <n; i++)
	  D[i] = INT_MAX;
  
	//initialize heap based on distance array of length n
	//heap has max size n (third argument)
	minHeapInit(D, n, n);

  D[s] = 0;
	//insert source vertex in heap
  minHeapInsert(s);
  while(!minHeapEmpty()){
		//extract min element from heap
		int i = heapExtractMin();
    for(int k=V[i]; k<V[i+1]; k++){
			int j = E[k];
			if(D[j] > D[i]+W[k]){
				D[j] = D[i] + W[k];
				if(isInHeap(j)){
					//update distance for vertex j and re-heapify
					minHeapChange(j);
				} else{
					//insert vertex j in heap
					minHeapInsert(j);
				}
			}
		}
  }
	for(int i=0;i<n;i++)
		printf("%d ",D[i]);
	printf("\n");
}

void readGraph(int *V, int *E, int *W, int n, int m){
	for(int i=0; i<n+1; i++)
		V[i] = 0;
	//create CSV arrays from edge list sorted by first vertex
	int vold=-1;
	for(int k=0; k<m; k++){
		int vi, vo, wt;
		if(scanf("%d %d %d", &vi, &vo, &wt)!= 3){
			fprintf(stderr, "input invalid\n");
			exit(1);
		}
		if(vi > n-1 || vo > n-1){
			fprintf(stderr, "vertex index too large\n");
			exit(1);
		}
		if(wt < 0){
			fprintf(stderr, "nonzero weights only\n");
			exit(1);
		}
		E[k] = vo;
		W[k] = wt;
		if(k == 0)
			V[vi] = 0;
		else if(vi != vold)
			V[vi] = k;
		vold = vi;
	}
	V[n] = m;
	//Find first out-edge
	int first = 0;
	for(int i = 1; i < n; i++)
		if(V[i] != 0){
			first = i-1;
			break;
		}
	//vertices with no out-edges
	for(int i = n-1; i > first; i--)
		if(V[i] ==0)
		 	V[i] = V[i+1];
}
