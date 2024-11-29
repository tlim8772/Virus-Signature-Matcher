
./gen_sig 1000 3000 3000 0.1 > sig.fasta
./gen_sample sig.fasta 980 20 1 2 10000 10000 10 30 0.1 > samp.fastq

./bench-a100 samp.fastq sig.fasta > ans
./matcher samp.fastq sig.fasta > out

diff ans out

if [ $? -ne 0 ]; then
    echo 'failed test'
else 
    rm sig.fasta samp.fastq ans out
fi 