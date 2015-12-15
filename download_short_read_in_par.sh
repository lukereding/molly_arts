# download sra files from ncbi short read archive
# first argument is the file name with the list of run ID's

echo $1

count=0
NR_CPUS=4

for i in $(cat "$1");
do
	printf $i
	prefetch $i &
	let count+=1
	echo $count
	[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
