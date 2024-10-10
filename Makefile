# Makefile

NVCC        := /usr/local/cuda/bin/nvcc
COMMON_SRCS := common.cc
HEADERS := kseq/kseq.h

CXXFLAGS := -std=c++20 -O3

OUTPUT_GEN  := gen_sample
OUTPUT_GEN_SIG  := gen_sig
OUTPUT_GPU  := matcher

.DEFAULT_TARGET: all

all: $(OUTPUT_GEN) $(OUTPUT_GEN_SIG) $(OUTPUT_GPU)

$(OUTPUT_GPU): kernel_skeleton.cu $(COMMON_SRCS)
	$(NVCC) $(CXXFLAGS) -lineinfo -dlto -arch=native -o $@ $^

$(OUTPUT_GEN): gen_sample.cc 
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OUTPUT_GEN_SIG): gen_sig.cc 
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OUTPUT_GEN) $(OUTPUT_GEN_SIG) $(OUTPUT_GPU) test_helper 

test_helper:
	$(NVCC) -arch=native test_helpers.cu -o test_helper

# for debugging matcher with cuda-gdb
debug_matcher:
	$(NVCC) -arch=native -g kernel_skeleton.cu common.cc -o debug_matcher

matcher:
	$(NVCC) -std=c++20 -O3 -lineinfo -dlto -arch=native kernel_skeleton.cu common.cc -o matcher