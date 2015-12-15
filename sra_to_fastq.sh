
count=0
NR_CPUS=4

for file in *.sra;
do
    printf $file
    fastq-dump --split-files $file &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done