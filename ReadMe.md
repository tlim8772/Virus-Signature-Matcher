if you compile on a100-80, and run on h100-96 (and vice versa), you will get errors.
testM: test match. test 3 is the worst case scenario of O(mn) time complexity. You can see the diff if you run on 1 thread vs many threads
testMS test match_score

i tested, and copying the data form host to device is not the bottleneck

if you use unified memory (i.e use cudaMallocManaged) you must use cudaDeviceSynchronise() after kernel call, this is because kernels run asynchronouslu
with respect to gpu

Don't use unified memory on SoC compute cluster, it is MUCH slower. Do cudaMemcpy as much as possible

compile matcher: `usr/local/bin/nvcc -arch=native -std=c++20 common.cc kernel_skeleton.cc -o matcher`

compile gen_sample: `g++ gen_sample.cc -o gen_sample`

compile gen_sig: `g++ gen_sig.cc -o gen_sig`

to run: `./matcher t1.fastq t1.fastq`