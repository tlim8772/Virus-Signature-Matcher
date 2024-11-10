#include "helpers.cu"
#include <chrono>
#include <iostream>
#include <string>

using namespace std;

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

        cudaMalloc(&phread_33, sizeof(char) * s.size());
        cudaMallocManaged(&res, sizeof(double));
        cudaMemcpy(phread_33, s.data(), sizeof(char) * s.size(), cudaMemcpyHostToDevice);

        //for (int i = 0; i < s.size(); i ++) phread_33[i] = s[i];
        *res = 0.0;

        int NUM_BLOCKS = (len_virus + BLOCK_SIZE) / BLOCK_SIZE;
        match_score<<<NUM_BLOCKS, BLOCK_SIZE>>>(phread_33, len_virus, idx, res);
        cudaDeviceSynchronize(); // important

        cout << *res << endl;
    }
}

void test_f() {
    int n;
    cin >> n;
    for (int i = 0; i < n; i ++) {
        string s, v;
        cin >> s >> v;

        char *sample;
        char *virus;
        int *res;

        cudaMalloc(&sample, sizeof(char) * s.size());
        cudaMalloc(&virus, sizeof(char) * v.size());
        cudaMallocManaged(&res, sizeof(int));

        *res = 999999;
        cudaMemcpy(sample, s.data(), s.size() * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(virus, v.data(), v.size() * sizeof(char), cudaMemcpyHostToDevice);

        //auto start = chrono::high_resolution_clock::now();
        
        f<<<1, 1>>>(sample, virus, s.size(), v.size(), res);
        cudaDeviceSynchronize();

        //auto end = chrono::high_resolution_clock::now();
        //auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        cout << (*res) << endl;
        //cout << "time: " << seconds.count() << endl;

    }
}
int main() {
    test_match();
    //test_match_score();
    //test_f();
}
