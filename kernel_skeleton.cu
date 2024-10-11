#include "kseq/kseq.h"
#include "common.h"
#include "helpers.cu"
#include "timer.cu"
#include <iostream>

#define MAX_LEN_SAMP 2200
#define MAX_LEN_SIG 1000
#define MAX 2200000
#define NUM_STREAMS 32

cudaStream_t streams[NUM_STREAMS];

char *samps[MAX_LEN_SAMP];
char *phread33[MAX_LEN_SAMP];
char *sigs[MAX_LEN_SIG];

// managed memory, so that cpu can access the results
int *matchIdxs;
double *scores;

int ROWS;
int COLS;

void allocMem(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures) {
    ROWS = samples.size(), COLS = signatures.size();
    
    // alloc mem for samples DNA seqs
    for (int i = 0; i < samples.size(); i ++) {
        cudaMalloc(&samps[i], sizeof(char) * samples[i].seq.size());
    }

    // alloc mem for phread33 
    for (int i = 0; i < samples.size(); i ++) {
        cudaMalloc(&phread33[i], sizeof(char) * samples[i].qual.size());
    }


    // alloc mem for signatures DNA seqs
    for (int i = 0; i < signatures.size(); i ++) {
        cudaMalloc(&sigs[i], sizeof(char) * signatures[i].seq.size());
    }

    // alloc matchIdxs
    cudaMallocManaged(&matchIdxs, sizeof(int) * MAX);
    for (int i = 0; i < MAX; i ++) matchIdxs[i] = 999999;

    // alloc score
    cudaMallocManaged(&scores, sizeof(double) * MAX);
    for (int i = 0; i < MAX; i ++) scores[i] = 0.0;
    
    // init streams
    for (int i = 0; i < NUM_STREAMS; i ++) {
        cudaStreamCreate(&streams[i]);
    }
}

void syncStreams() {
    for (int i = 0; i < NUM_STREAMS; i ++) cudaStreamSynchronize(streams[i]);
}


void runMatcher(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures, std::vector<MatchResult>& matches) {
   
    allocMem(samples, signatures);
    
   
    for (int i = 0; i < ROWS; i ++) {
        //int sid = i & (NUM_STREAMS - 1);
        cudaMemcpy(samps[i], samples[i].seq.data(), samples[i].seq.size() * sizeof(char), cudaMemcpyHostToDevice);
    }
    for (int j = 0; j < COLS; j ++) {
        //int sid = j & (NUM_STREAMS - 1);
        cudaMemcpy(sigs[j], signatures[j].seq.data(), signatures[j].seq.size() * sizeof(char), cudaMemcpyHostToDevice);
    }
    
    //syncStreams();
    
    for (int i = 0; i < ROWS; i ++) {
        for (int j = 0; j < COLS; j ++) {
            int idx = i * ROWS + j;
            int sid = idx & (NUM_STREAMS - 1); // take idx % 32
            int sampLen = samples[i].seq.size();
            int sigLen = signatures[j].seq.size();
            int numBlocks = (sampLen + BLOCK_SIZE) / BLOCK_SIZE;
            
            match<<<numBlocks, BLOCK_SIZE, 0, streams[sid]>>>(samps[i], sigs[j], sampLen, sigLen, &matchIdxs[idx]);
        }
    }

    syncStreams();
    cudaDeviceSynchronize();

    bool copied[ROWS] = {false};
    for (int i = 0; i < ROWS; i ++) {
        for (int j = 0; j < COLS; j ++) {
            int idx = i * ROWS + j;
            //int sid = idx & (NUM_STREAMS - 1); // take idx % 32

            // if pattern is found
            if (matchIdxs[idx] < INVALID) {
                int sigLen = signatures[j].seq.size();
                int numBlocks = (sigLen + BLOCK_SIZE) / BLOCK_SIZE;

                if (!copied[i]) {
                    copied[i] = true;
                    cudaMemcpy(phread33[i], samples[i].qual.data(), samples[i].qual.size() * sizeof(char), cudaMemcpyHostToDevice);
                }
                
                match_score<<<numBlocks, BLOCK_SIZE>>>(phread33[i], sigLen, matchIdxs[idx], &scores[idx]);
            }
        }
    }

    //syncStreams();
    cudaDeviceSynchronize();
    for (int i = 0; i < ROWS; i ++) {
        for (int j = 0; j < COLS; j ++) {
            int idx = i * ROWS + j;
            if (matchIdxs[idx] < INVALID) matches.push_back({samples[i].name, signatures[j].name, scores[idx]});
        }
    }    
}
