#include <bits/stdc++.h>
#include <string.h>
#include <cuda_runtime.h>
#include <chrono>

#define BLOCK_SIZE 128
#define INVALID 999999

using namespace std;


// len is the actual length of the string, and does not include the null character at the end
// for string sample and virus, find the leftmost matching index
// each thread check for suffix of sample starting at tid 
__global__ void match(char* sample, char* virus, int len_sample, int len_virus, int* res) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ int matched[BLOCK_SIZE];
    
    // matched[threadIdx.x] is the idx of sample if pattern matters at idx
    if (tid >= len_sample || len_virus > len_sample - tid) {
        matched[threadIdx.x] = INVALID;
    } else {
        int i = 0;
        for (; i < len_virus; i ++) {
            if (virus[i] != 'N' && sample[i + tid] != 'N' && virus[i] != sample[i + tid]) break;
        }
        matched[threadIdx.x] = (i == len_virus) ? tid : INVALID;
    }
    
    
    // get the min of all idx of sample that has a match with virus
    for (int i = 2; i <= BLOCK_SIZE; i <<= 1) {
        __syncthreads();
        // if threadIdx.x is a multiple of i
        if (!(threadIdx.x & (i - 1))) {
            matched[threadIdx.x] = min(matched[threadIdx.x], matched[threadIdx.x + (i >> 1)]);
        }
        
    }
    if (!threadIdx.x) atomicMin(res, matched[0]);

}

// get the match score
// idx (of sample) is where the pattern is found to match
__global__ void match_score(char* phread_33, int len_virus, int idx, double* res) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ double accu[BLOCK_SIZE];

    accu[threadIdx.x] = (tid >= len_virus) ? 0.0 : (phread_33[idx + tid] - 33) / (double) len_virus;
    

    for (int i = 2; i <= BLOCK_SIZE; i <<= 1) {
        __syncthreads();
        // if threadIdx.x is a multiple of i
        if (!(threadIdx.x & (i - 1))) {
            accu[threadIdx.x] += accu[threadIdx.x + (i >> 1)];
        }
       
    }
    
    if (!threadIdx.x) atomicAdd(res, accu[0]);
}
