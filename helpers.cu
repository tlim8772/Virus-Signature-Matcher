#include <bits/stdc++.h>
#include <string.h>
#include <cuda_runtime.h>
#include <chrono>

#define BLOCK_SIZE 128
#define INVALID 999999



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
__global__ void match_score(char* phread33, int len_virus, int idx, double* res) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ double accu[BLOCK_SIZE];
    
    accu[threadIdx.x] = (tid >= len_virus) ? 0.0 : (phread33[idx + tid] - 33) / (double) len_virus;
    
    for (int i = 2; i <= BLOCK_SIZE; i <<= 1) {
        __syncthreads();
        // if threadIdx.x is a multiple of i
        if (!(threadIdx.x & (i - 1))) {
            accu[threadIdx.x] += accu[threadIdx.x + (i >> 1)];
        }
       
    }
    
    if (!threadIdx.x) atomicAdd(res, accu[0]);
}

// to test time taken for 1 thread only to do naive matching
__global__ void f(char* samp, char* sig, int sampLen, int sigLen, int* res) {
    
    for (int i = 0; i <= sampLen - sigLen; i ++) {
        int j = 0;
        for (; j < sigLen; j ++) {
            if (samp[i + j] != 'N' && sig[j] != 'N' && samp[i + j] != sig[j]) break;
        }
        if (j == sigLen) {
            *res = i;
            break;
        }
    }
}

__global__ void matcherKernel(char **samps, char **sigs, char **phread33s, int *sampLens, int *sigLens, int ROWS, int COLS, double *results) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ROWS * COLS) return;

    int r = tid / COLS;
    int c = tid % COLS;
    char *samp = samps[r];
    char *sig = sigs[c];
    int sampLen = sampLens[r];
    int sigLen = sigLens[c];

    int lim = sampLen - sigLen;
    for (int i = 0; i <= lim; i ++) {
        int j = 0;
        for(; j < sigLen; j ++) {
            if (samp[i + j] != 'N' && sig[j] != 'N' && samp[i + j] != sig[j]) break;
        }
        
        // if we are successful, calculate the score
        if (j == sigLen) {
            char *phread33 = phread33s[r];
            int score = 0;
            for (int k = 0; k < sigLen; k ++) score += (phread33[i + k] - 33);
            results[tid] = score / (double) sigLen;
            return;
        }
    }
    
}


// combine both methods
// each pair gets a block, with has 128 threads to split up the work, so we have a 2d grid of blocks
__global__ void combineBoth(char **samps, char **sigs, char **phread33s, int *sampLens, int *sigLens, int ROWS, int COLS, double *score) {
    int r = blockIdx.x;
    int c = blockIdx.y;
    

    char *samp = samps[r];
    char *sig = sigs[c];
    char *phread33 = phread33s[r];
    int sampLen = sampLens[r];
    int sigLen = sigLens[c];
    int lim = sampLen - sigLen;

    __shared__ int store[BLOCK_SIZE];
    
    store[threadIdx.x] = 999999;

    __syncthreads();

    for (int i = 0; i <= lim; i += BLOCK_SIZE) {
        if (i + threadIdx.x > lim) continue;

        int idx = i + threadIdx.x;
        int j = 0;
        for (;j < sigLen; j ++) {
            if (samp[idx  + j] != 'N' && sig[j] != 'N' && samp[idx + j] != sig[j]) break;
        }

        if (j == sigLen) store[threadIdx.x] = min(store[threadIdx.x], idx);
    }

    // coalse the result
    for (int i = 2; i <= BLOCK_SIZE; i <<= 1) {
        __syncthreads();
        
        if (!(threadIdx.x & (i - 1))) {
            store[threadIdx.x] = min(store[threadIdx.x], store[threadIdx.x + (i >> 1)]);
        }
       
    }

    __syncthreads();
    if (store[0] == 999999) return;
    
   int matchIdx = store[0];
    __syncthreads();

    store[threadIdx.x] = 0;
    __syncthreads();

    for (int i = 0; i < sigLen; i += BLOCK_SIZE) {
        store[threadIdx.x] += (i + threadIdx.x >= sigLen) ? 0 : (phread33[matchIdx + i + threadIdx.x] - 33);
    }

    for (int i = 2; i <= BLOCK_SIZE; i <<= 1) {
        __syncthreads();
        
        if (!(threadIdx.x & (i - 1))) {
            store[threadIdx.x] += store[threadIdx.x + (i >> 1)];
        }
       
    }
    

    if (threadIdx.x == 0) score[r * COLS + c] = store[0] / (double) sigLen;
}