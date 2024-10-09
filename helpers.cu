#include <bits/stdc++.h>
#include <string.h>
#include <cuda_runtime.h>
#include <chrono>

#define BLOCK_SIZE 256
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
        // if threadIdx.x is a multiple of i
        if (!(threadIdx.x & (i - 1))) {
            matched[threadIdx.x] = min(matched[threadIdx.x], matched[threadIdx.x + (i >> 1)]);
        }
        __syncthreads();
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
        // if threadIdx.x is a multiple of i
        if (!(threadIdx.x & (i - 1))) {
            accu[threadIdx.x] += accu[threadIdx.x + (i >> 1)];
        }
        __syncthreads();
    }
    
    if (!threadIdx.x) atomicAdd(res, accu[0]);
}


void test_match() {
    int n;
    cin >> n;

    for (int i = 0; i < n; i ++) {
        string s, v;
        cin >> s >> v;

        char *sample;
        char *virus;
        int *res;

        cudaMallocManaged(&sample, sizeof(char) * s.size());
        cudaMallocManaged(&virus, sizeof(char) * v.size());
        cudaMallocManaged(&res, sizeof(int));

        for (int i = 0; i < s.size(); i ++) sample[i] = s[i];
        for (int i = 0; i < v.size(); i ++) virus[i] = v[i];
        *res = INVALID;

        int NUM_BLOCKS = (s.size() + BLOCK_SIZE) / BLOCK_SIZE;

        auto start = chrono::high_resolution_clock::now();
        match<<<NUM_BLOCKS, BLOCK_SIZE>>>(sample, virus, s.size(), v.size(), res);
        cudaDeviceSynchronize(); // important

        
        auto end = chrono::high_resolution_clock::now();

        
        auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        cout << *res << endl;
        cout << "Time taken (ms): " << seconds.count() << endl;
    }
}

void test_match_score() {
    int n;
    cin >> n;
    for (int i = 0; i < n; i ++) {
        string s;
        int len_virus, idx;
        cin >> s >> len_virus >> idx;

        char *phread_33;
        double *res;

        cudaMallocManaged(&phread_33, sizeof(char) * s.size());
        cudaMallocManaged(&res, sizeof(double));

        for (int i = 0; i < s.size(); i ++) phread_33[i] = s[i];
        *res = 0.0;

        int NUM_BLOCKS = (len_virus + BLOCK_SIZE) / BLOCK_SIZE;
        match_score<<<NUM_BLOCKS, BLOCK_SIZE>>>(phread_33, len_virus, idx, res);
        cudaDeviceSynchronize(); // important

        cout << *res << endl;
    }
}

int main() {
    //test_match();
    test_match_score();
}