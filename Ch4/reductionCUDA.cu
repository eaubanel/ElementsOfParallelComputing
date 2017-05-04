/* Copyright 2017 Eric Aubanel
 * This file contains code implementing Algorithm 4.14 from
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
 * Implementation of reduction of n floats using CUDA
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int isPowerOf2(int n);

//reduce n floats in array a to a partial sum for each block,
//stored in array c. Block size must be power of 2
__global__ void reductionGPU(float *a, float *c, int n){
	//size of b indicated by kernel call in main
  extern __shared__ float b[];
	int gsize = blockDim.x; //block size
	int nt = gsize * gridDim.x; //total number of threads
	int gid = blockIdx.x; //block id
	int tid = threadIdx.x; //local thread id
	int id = gid*gsize + tid; //global thread id
	//evaluate as float to avoid overflow
	int istart = (float)id*n/nt;
	int iend = (float)(id+1)*n/nt - 1;
	
	float psum = 0.0;
	for(int i=istart; i<=iend; i++)
		psum += a[i];
	b[tid] = psum;

  __syncthreads();
	for(int j=gsize>>1; j>=1; j >>= 1){
		if(tid<j)
			b[tid] += b[tid+j];
		__syncthreads();
	}
	c[gid] = b[0];
}

int main(int argc, char **argv){
	float *a_h; //array to be reduced on host (CPU)
	float *c_h; //array of partial sums on host
	float *a_d; //array to be reduced on device (GPU)
	float *c_d; //array of partial sums on device
	cudaError_t error1, error2;
	struct timespec tstart, tend;
	float time;

	if(argc < 4){
		fprintf(stderr,"usage: %s n blockSize numBlocks\n", argv[0]);
		return 1;
	}
	int n = strtol(argv[1], NULL, 10);
	int blockSize = strtol(argv[2], NULL, 10); //size of thread block on device
	int numBlocks = strtol(argv[3], NULL, 10); //number of blocks on device
	if(!isPowerOf2(blockSize)){
		fprintf(stderr,"blockSize must be power of 2\n");
		return 1;
	}

	//memory allocation on host and device
	a_h = (float *)malloc(n*sizeof(float));
	c_h = (float *)malloc(numBlocks*sizeof(float));
	if(a_h == NULL || c_h == NULL){
		fprintf(stderr,"couldn't allocate memory on host\n");
		return 1;
	}
	error1 = cudaMalloc((void **)&a_d, n*sizeof(float));
	error2 = cudaMalloc((void **)&c_d, numBlocks*sizeof(float));
	if(error1 != cudaSuccess || error2 != cudaSuccess){
		fprintf(stderr,"couldn't allocate memory on device\n");
		return 1;
	}

	for(int i=0; i<n; i++)
		a_h[i] = rand()%100;
	//sequential reduction for verification and timing
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	float sum = 0.0;
	for(int i=0; i<n; i++)
		sum += a_h[i];
	clock_gettime(CLOCK_MONOTONIC, &tend);
	time = (tend.tv_sec-tstart.tv_sec) + (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("CPU reduction time in s: %f\n", time);

	//timing won't include transfer of array to device
	cudaMemcpy(a_d, a_h, n*sizeof(float), cudaMemcpyHostToDevice);
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	reductionGPU <<<numBlocks, blockSize, blockSize*sizeof(float)>>> (a_d, c_d, n);
	error1 = cudaThreadSynchronize();// wait until GPU kernel finished
	if(error1 != cudaSuccess){
		fprintf(stderr,"error executing kernel: %s\n", cudaGetErrorString(error1));
		return 1;
	}
	cudaMemcpy(c_h, c_d, numBlocks*sizeof(float), cudaMemcpyDeviceToHost);
	float dsum =0.0;
	for(int i=0; i<numBlocks; i++)
		dsum += c_h[i];
	clock_gettime(CLOCK_MONOTONIC, &tend);
	time = (tend.tv_sec-tstart.tv_sec) + (tend.tv_nsec-tstart.tv_nsec)*1.0e-9;
	printf("GPUS time in s: %f\n", time);

	//not necessarily the same because of differences in roundoff error
	printf("relative difference between sequential and parallel sums: %g\n",
			fabs(dsum-sum)/sum);
	return 0;
}

int isPowerOf2(int n){
	while(n){
		if(n & 1)
			break;
		n >>= 1;
	}
	return (1 == n? 1:0);
}
