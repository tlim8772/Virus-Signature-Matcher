#include "kseq/kseq.h"
#include "common.h"
#include "helpers.cu"
#include "timer.cu"
#include <iostream>

char **samps, **dsamps;
char **phread33, **dphread33;
char **sigs, **dsigs;

double *scores, *dscores;
int *sampLens, *dsampLens;
int *sigLens, *dsigLens;

int MAX;
int ROWS;
int COLS;

void allocMem(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures) {
    samps = (char**) malloc(sizeof(char*) * ROWS);
    phread33 = (char**) malloc(sizeof(char*) * ROWS);
    sigs = (char**) malloc(sizeof(char*) * COLS);
    
    cudaMalloc((void**)&dsamps, sizeof(char*) * ROWS);
    cudaMalloc((void**)&dphread33, sizeof(char*) * ROWS);
    cudaMalloc((void**)&dsigs, sizeof(char*) * COLS);

    scores = (double*) malloc(sizeof(double) * MAX);
    sampLens = (int*) malloc(sizeof(int) * ROWS);
    sigLens = (int*) malloc(sizeof(int) * COLS);

    cudaMalloc(&dscores, sizeof(double) * MAX);
    cudaMalloc(&dsampLens, sizeof(int) * ROWS);
    cudaMalloc(&dsigLens, sizeof(int) * COLS);
    
    
    for (int i = 0; i < ROWS; i ++) {
        cudaMalloc(&samps[i], sizeof(char) * samples[i].seq.size());
    }

    
    for (int i = 0; i < ROWS; i ++) {
        cudaMalloc(&phread33[i], sizeof(char) * samples[i].qual.size());
    }

    for (int i = 0; i < COLS; i ++) {
        cudaMalloc(&sigs[i], sizeof(char) * signatures[i].seq.size());
    }

   

    for (int i = 0; i < MAX; i ++) scores[i] = -999999.0;

    for (int i = 0; i < ROWS; i ++) sampLens[i] = samples[i].seq.size();

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
    
    cudaMemcpy(dsamps, samps, sizeof(char*) * ROWS, cudaMemcpyHostToDevice);
    cudaMemcpy(dphread33, phread33, sizeof(char*) * ROWS, cudaMemcpyHostToDevice);
    cudaMemcpy(dsigs, sigs, sizeof(char*) * COLS, cudaMemcpyHostToDevice);
    
    cudaMemcpy(dscores, scores, sizeof(double) * MAX, cudaMemcpyHostToDevice);
    cudaMemcpy(dsampLens, sampLens, sizeof(int) * ROWS, cudaMemcpyHostToDevice);
    cudaMemcpy(dsigLens, sigLens, sizeof(int) * COLS, cudaMemcpyHostToDevice);
    

    int NUM_BLKS = (MAX + BLOCK_SIZE) / BLOCK_SIZE;
    matcherKernel<<<NUM_BLKS, BLOCK_SIZE>>>(dsamps, dsigs, dphread33, dsampLens, dsigLens, ROWS, COLS, dscores);
    cudaMemcpy(scores, dscores, sizeof(double) * MAX, cudaMemcpyDeviceToHost);
    //cudaDeviceSynchronize();
    
    for (int i = 0; i < ROWS; i ++) {
        for (int j = 0; j < COLS; j ++) {
            int idx = i * COLS + j;
            if (scores[idx] > -0.000001) {
                matches.push_back({samples[i].name, signatures[j].name, scores[idx]});
            }
        }
    }
}
