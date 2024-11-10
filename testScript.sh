for i in {1..9}; do
    ./matcher tests/t$i.fastq tests/t$i.fasta > t$i.ans
    diff t$i.ans tests/t$i.out
    
    if [ $? -ne 0 ]; then
        echo "Failed on test case $i"
        break
    else 
        rm t$i.ans
    fi

done