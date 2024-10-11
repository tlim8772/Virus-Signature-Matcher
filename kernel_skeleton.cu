#include "kseq/kseq.h"
#include "common.h"
#include "helpers.cu"
#include "timer.cu"
#include <iostream>

char **samps;
char **phread33;
char **sigs;

double *scores;
int *sampLens;
int *sigLens;

int MAX;
int ROWS;
int COLS;

void allocMem(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures) {
    cudaMallocManaged((void**) &samps, sizeof(char*) * ROWS);
    cudaMallocManaged((void**) &phread33, sizeof(char*) * ROWS);
    cudaMallocManaged((void**) &sigs, sizeof(char*) * COLS);
    
    // alloc mem for samples DNA seqs
    for (int i = 0; i < ROWS; i ++) {
        cudaMalloc(&samps[i], sizeof(char) * samples[i].seq.size());
    }

    // alloc mem for phread33 
    for (int i = 0; i < ROWS; i ++) {
        cudaMalloc(&phread33[i], sizeof(char) * samples[i].qual.size());
    }


    // alloc mem for signatures DNA seqs
    for (int i = 0; i < COLS; i ++) {
        cudaMalloc(&sigs[i], sizeof(char) * signatures[i].seq.size());
    }

    // alloc score
    cudaMallocManaged(&scores, sizeof(double) * MAX);
    for (int i = 0; i < MAX; i ++) scores[i] = -999999.0;

    cudaMallocManaged(&sampLens, sizeof(int) * ROWS);
    for (int i = 0; i < ROWS; i ++) sampLens[i] = samples[i].seq.size();

    cudaMallocManaged(&sigLens, sizeof(int) * COLS);
    for (int i = 0; i < COLS; i ++) sigLens[i] = signatures[i].seq.size();
    
}



void runMatcher(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures, std::vector<MatchResult>& matches) {
    ROWS = samples.size(), COLS = signatures.size();
    MAX = ROWS * COLS;
    
    allocMem(samples, signatures);
    
    for (int i = 0; i < ROWS; i ++) {
        cudaMemcpy(samps[i], samples[i].seq.data(), sizeof(char) * samples[i].seq.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(phread33[i], samples[i].qual.data(), sizeof(char) * samples[i].qual.size(), cudaMemcpyHostToDevice);
    }
    for (int i = 0; i < COLS; i ++) {
        cudaMemcpy(sigs[i], signatures[i].seq.data(), sizeof(char) * signatures[i].seq.size(), cudaMemcpyHostToDevice);
    }

    int NUM_BLKS = (MAX + BLOCK_SIZE) / BLOCK_SIZE;
    matcherKernel<<<NUM_BLKS, BLOCK_SIZE>>>(samps, sigs, phread33, sampLens, sigLens, ROWS, COLS, scores);
    cudaDeviceSynchronize();
    
    for (int i = 0; i < ROWS; i ++) {
        for (int j = 0; j < COLS; j ++) {
            int idx = i * COLS + j;
            if (scores[idx] > -0.000001) {
                matches.push_back({samples[i].name, signatures[j].name, scores[idx]});
            }
        }
    }
}
