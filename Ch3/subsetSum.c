/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm  3.6 from
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
 * Implementation of T/F subset sum problem using dynamic programming
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int main(int argc, char **argv){
  int R; // max magnitude of elements in set
  int n; // number of points in set
	int S; // target sum = nR/4
	int *s; // set of points (starts at index 1)
	char *F; // Dynamic programming table

	if(argc < 3){
		fprintf(stderr,"usage: %s R n [seed]\n", argv[0]);

		return 1;
	}
	R = strtol(argv[1], NULL, 10);
	n = strtol(argv[2], NULL, 10);
  S = n*R/4;
  int m = S+1;
  s = malloc((n+1)*sizeof(int));
  F = calloc((n+1)*m,sizeof(char));
	if(!s || !F){
		fprintf(stderr,"couldn't allocate memory\n");
		return 1;
	}
	if(4 == argc) 
		srand(strtol(argv[3], NULL, 10));
	else
		srand(time(NULL));
  for(int i = 1; i<=n;i++){
		s[i] = rand()%R;
		printf("%d ", s[i]);
		F[i*m] = 1;
  }
	printf("\n sum = %d\n", S);
  F[m+s[1]] = 1;
  for(int i=2;i<=n;i++){
    for(int j=1; j<s[i];j++){
      F[i*m+j] = F[(i-1)*m+j];
    }
    for(int j=s[i];j<=S;j++){
      F[i*m+j] = F[(i-1)*m+j] || F[(i-1)*m+j-s[i]];
    }
  }
  printf("%s\n",F[n*m+S]?"true":"false");
  return 0;
}
