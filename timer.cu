#include <string>

struct Timer {
    cudaEvent_t begin, end;

    Timer() {
        cudaEventCreate(&begin);
        cudaEventCreate(&end);
    }

    ~Timer() {
        cudaEventDestroy(begin);
        cudaEventDestroy(end);
    }


    void start() {
        cudaEventRecord(begin);
    }

    void stop() {
        cudaEventRecord(end);
    }

    std::string time() {
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, begin, end);
        return std::to_string(milliseconds / 1000);
    }

};