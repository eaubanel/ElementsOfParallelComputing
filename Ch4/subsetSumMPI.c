/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithms 4.19 & 4.20 from
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
 * Parallel implementation of T/F subset sum problem using dynamic programming
 * and MPI.
 * Result validated by comparing to sequential solution.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

#define MAX(x, y) ((x)>(y)?(x):(y))
#define FIND_ID(j, p, n) (((p)*((j)+1)-1)/(n))
#define CEIL(a, b) (((a)+(b)-1)/(b))

//Sequential solution
char *solveSequential(int S, int *s, int n);
//Gather result array from all processes
char *gatherParallel(char *F, int S, int n, int nb, int p, int id);
//Solve row i of F with array L for values obtained from other processes
void solveRow(char *F, char*L, int *s, int nb, int myFirst, int i, int id);

int main(int argc, char **argv){
  int R; // max magnitude of elements in set
  int n; // number of points in set
	int S; // target sum = nR/4
	int *s; // set of points (starts at index 1)
	char *Fs; // Dynamic programming table (sequential)
	char *F; // Dynamic programming table (distributed)
	char *Fp;// Dynamic programming table (gathered at process 0)
	char *L; //buffer for receives
	int id; //my id
	int p; //number of processes
	double timer;

	MPI_Status status;
  MPI_Request req1, req2;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if(argc < 3){
		if(!id) fprintf(stderr,"usage: %s R n [seed]\n", argv[0]);
		MPI_Finalize();
		return 1;
	}
	R = strtol(argv[1], NULL, 10);
	n = strtol(argv[2], NULL, 10);
  S = n*R/4;
	int myFirst = id*S/p + 1;
	int myLast = (id+1)*S/p;
	if(!id)
		myFirst = 0;
	int nb = myLast - myFirst + 1;
	int nL = CEIL(S, p);
  s = malloc((n+1)*sizeof(int));
  F = calloc((n+1)*nb,sizeof(char));
	L = malloc((S/p)*sizeof(char));
	if(!s || !F || !L){
		if(!id) fprintf(stderr,"couldn't allocate memory\n");
		MPI_Finalize();
		return 1;
	}
	if(!id)
		F[nb] = 1;

	if(!id){
		if(4 == argc) 
			srand(strtol(argv[3], NULL, 10));
		else
			srand(time(NULL));
  	for(int i = 1; i<=n;i++){
			s[i] = rand()%R;
			printf("%d ", s[i]);
		}
		printf("\n sum = %d\n", S);
  }
	MPI_Bcast(s, n+1, MPI_INT, 0, MPI_COMM_WORLD);
		
	if(id == FIND_ID(s[1]-1, p, S))
  	F[nb+s[1]-myFirst] = 1;

	MPI_Barrier(MPI_COMM_WORLD);
	timer = -MPI_Wtime();

  for(int i=2; i<=n; i++){
		int id1 = FIND_ID(myFirst+s[i]-1, p, S);
		int id2 = FIND_ID(myLast+s[i]-1, p, S);
		int myLocalBegin;
		if(id1 < p){
			if(id1 == id2){
				myLocalBegin = MAX(0,id1*S/p+1-s[i]-myFirst);
				MPI_Isend(F+(i-1)*nb+myLocalBegin, nb-myLocalBegin, MPI_CHAR, id1, 0, 
						MPI_COMM_WORLD, &req1);
			}else{
				int destBegin = myFirst + s[i];
				int destLast = (id1+1)*S/p;
				int nb1 = destLast - destBegin + 1;
				if(id1 > id) //bug at line 22 in Algorithm 4.19
					MPI_Isend(F+(i-1)*nb, nb1, MPI_CHAR, id1, 0,
							MPI_COMM_WORLD, &req1);
				if(id2 < p){
					MPI_Isend(F+(i-1)*nb+nb1, nb-nb1, MPI_CHAR, id2, 0,
							MPI_COMM_WORLD, &req2);
				}
			}
		}
		if(id && (myLast - s[i] >= 0)){ //bug at line 26 of Algorithm 4.19
			id1 = FIND_ID(myFirst-s[i]-1, p, S);
			if(myFirst - s[i] < 0)
				myLocalBegin = s[i] - myFirst;
			else
				myLocalBegin = 0;
			MPI_Recv(L+myLocalBegin, nL-myLocalBegin, MPI_CHAR, id1, 0,
					MPI_COMM_WORLD, &status);
			int nS;
			MPI_Get_count(&status, MPI_CHAR, &nS);
			id2 = FIND_ID(myLast - s[i] - 1, p, S);
			if(id1 != id2 && id2 < id){ //bug at line 36 of Algorithm 4.19
				MPI_Recv(L+nS, nL-nS, MPI_CHAR, id2, 0,
					MPI_COMM_WORLD, &status);
			}
		}
		solveRow(F, L, s, nb, myFirst, i, id);
	}
	timer += MPI_Wtime();
	double ptime;
	MPI_Reduce(&timer, &ptime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(!id)
		printf("time in seconds: %f\n", ptime);
	if(id == p-1)
		printf("%s\n",F[n*nb+nb-1]?"true":"false");

	//sequential solution
	if(!id){
		Fs = solveSequential(S, s, n);
		if(!Fs){
			fprintf(stderr,"couldn't allocate memory for sequential solution\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	//gather parallel solution and verify it
	Fp = gatherParallel(F, S, n, nb, p, id);
	if(!id){
		int m = S+1;
		int passed = 1;
		for(int i=2;i<=n;i++){
    	for(int j=1; j<=S;j++){
				if(Fs[i*m+j] != Fp[i*m+j]){
					printf("i=%d, j=%d, Fs=%d, Fp=%d\n", i, j, Fs[i*m+j], Fp[i*m+j]);
					passed = 0;
				}
			}
		}
		if(passed)
			printf("result verified\n");		
	}
	MPI_Finalize();
  return 0;
}

char *solveSequential(int S, int *s, int n){
	int m = S+1;
	char *Fs = calloc((n+1)*m,sizeof(char));
	if(Fs){
		for(int i = 1; i<=n;i++)
			Fs[i*m] = 1;
		Fs[m+s[1]] = 1;
		for(int i=2;i<=n;i++){
    	for(int j=1; j<s[i];j++){
      	Fs[i*m+j] = Fs[(i-1)*m+j];
    	}
    	for(int j=s[i];j<=S;j++){
      	Fs[i*m+j] = Fs[(i-1)*m+j] || Fs[(i-1)*m+j-s[i]];
   		}
  	}
	}
	return Fs;
}

char *gatherParallel(char *F, int S, int n, int nb, int p, int id){
	char *Fp;
	if(!id){
		Fp = malloc((n+1)*(S+1)*sizeof(char));
		if(!Fp){
			fprintf(stderr,"couldn't allocate memory for validation\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	//Create arrays needed for Gatherv;
	//disp and recvCount arrays only used by process 0
	int disp[p];
	int recvCount[p];
	if(!id){
		recvCount[0] = S/p + 1;
		disp[0] = 0;
		for(int i=1; i<p; i++){
			recvCount[i] = (i+1)*S/p - i*S/p;
			disp[i] = disp[i-1]+ recvCount[i-1];
		}
	}
	//gather parallel result, row-by-row
	for(int i=2; i<=n; i++){
		MPI_Gatherv(F+i*nb, nb, MPI_CHAR, Fp+i*(S+1), recvCount, disp, MPI_CHAR, 0,
				MPI_COMM_WORLD);
	}
	return Fp;
}


void solveRow(char *F, char*L, int *s, int nb, int myFirst, int i, int id){
	int jstart = 0;
	if(!id){
		F[i*nb] = 1;
		jstart = 1;
	} 
	for(int j=jstart; j<nb; j++){
  	F[i*nb+j] = F[(i-1)*nb+j];
		int offset = j - s[i];
		if(offset >= 0)
    	F[i*nb+j] = F[(i-1)*nb+j] || F[(i-1)*nb+offset];
		else if(offset + myFirst >= 0)
			F[i*nb+j] = F[(i-1)*nb+j] || L[j];
  }
}
