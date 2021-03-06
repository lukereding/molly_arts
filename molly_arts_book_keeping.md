# playing around with sailfin RNA-seq data from Fraser et al. 2014
## getting the data
The data are stored on NCBI's short read archive, which is a terrible site. You cannot press a button or copy a link to download the raw data. I had to follow instructions [here](https://www.biostarhandbook.com/unit/sra/short-read-archive.html#sra). I first installed _sratoolkit_ on my Mac using Homebrew like `brew install sratoolkit`. I could then use the `prefetch` command to download the reads corresponding to each individual. I put the names of all the runs in a text document called runs, then I wrote a script, `download_short_read_in_par.sh`, to download the data in parrallel (this script and others are in this repository). Once the data are downloaded, they are in `.sra` format and not `.fastq` like you need them to be. So I wrote a second script, `sra_to_fastq.sh`, to convert each `.sra` file to `.fastq`, which takes a bit of time. Finally, I transferred all 11 fastq files to lonestar using `scp /path/to/fastq/files/*.fastq reding@lonestar.tacc.utexas.edu`.         

The reads are 100 bp single end reads. `head -4 SRR1161450_1.fastq` gives:

>@SRR1161450.1 HWI-ST619:197:D0VKFACXX:6:1101:1073:1995 length=100 CAGNTTTAGTCCAAAGTTTCTATATACAGTCAGAGATGAAACAGTTCTGGGCTTGGCCAAGCTGAAAAGAGGCTTCAGCTCCAGCTGAGTTCATCATTTN +SRR1161450.1 HWI-ST619:197:D0VKFACXX:6:1101:1073:1995 length=100 B@C#4ADDHHHHHIJJHIJIJFIIJJJIJIJJJJJIJJHIJJJJGIJJJJJJJJJJJJJJJJJJHGIJJIJJHFHHFFFFFFEECEECB;CAACDDDEC#

**Note**: It looks like there is a way to do all this just on TACC (not using your computer to download everything), although it would require building the software on TACC, which doesn't seem like it would be fun. There are brief instructions on how to install custom software to TACC under the "installing custom software" header on [this](https://wikis.utexas.edu/display/CoreNGSTools/Running+batch+jobs+at+TACC) page.

--------------------------------------------------------------------------------

## quality filtering
In the paper, they don't say if or how they did quality filtering. Totally arbitrarily, I decided to filter the reads so that at least 80% of each read has a quality score of 30 or above.

`pwd`<br>/scratch/02535/reding/molly_arts     

#### creating a job file
`touch filter_job`<br>`for i in *.fastq; do echo "fastq_quality_filter -q 20 -p 80 -i $i -Q 33 -o $i.filtered" >> filter_job; done`<br>`cat filter_job`       

>fastq_quality_filter -q 20 -p 80 -i SRR1161450_1.fastq -Q 33 -o SRR1161450_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1161451_1.fastq -Q 33 -o SRR1161451_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1165201_1.fastq -Q 33 -o SRR1165201_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1165203_1.fastq -Q 33 -o SRR1165203_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166366_1.fastq -Q 33 -o SRR1166366_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166367_1.fastq -Q 33 -o SRR1166367_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166368_1.fastq -Q 33 -o SRR1166368_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166369_1.fastq -Q 33 -o SRR1166369_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166370_1.fastq -Q 33 -o SRR1166370_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166371_1.fastq -Q 33 -o SRR1166371_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166372_1.fastq -Q 33 -o SRR1166372_1.fastq.filtered

That should work fine as a job file.

`module load fastx_toolkit` `module load launcher`       

#### submiting the job
From here on out, I'm following instructions [here](https://wikis.utexas.edu/display/bioiteam/Submitting+Jobs+to+Lonestar) to submit jobs to Lonestar.      

`cp $TACC_LAUNCHER_DIR/launcher.sge ./`      

I changed the laucher script so that the parameters looked like this:          <------ Setup Parameters ------>

$ -N Parametric $ -pe 12way 12 $ -q normal $ -o Parametric.o$JOB_ID $ -l h_rt=6:00:00 $ -V $ -M lukereding@utexas.edu. $ -m be $ -cwd

Actually, fuck that, let's use the python script to create our script instead:

`module load python`<br>`launcher_creator.py -n filter_reads -a Sailfin_RNASeq -q normal -t 6:00:00 -e lukereding@utexas.edu -j filter_job`<br>`qsub launcher.sge`       

To check on the status of your job, use `qstat`. When it's done running, you should have 11 fastq.filtered files: `ls *.fastq.filtered | wc -l`.      

### looking at the results of filtering
To see how many reads we weeded out in each of the samples, first write one line for each of the 11 original files:
`for i in *.fastq; do echo "orginally: grep ^@ $i | wc -l , filtered: grep ^@ $i".filtered" | wc -l >> filtering_results" >> compare_filtered; done`.
Note that we want bash to evaluate the grep expression in each line; we can tell it to that using backticks like `evaluate this expression`. This series of sed substituitons puts that backticks in place: ``cat compare_filtered  | sed 's,grep,`grep,g' | sed 's,wc -l, wc -l`,g' > filt_vs_original``.
Finally, we have to wrap each line in a call to `echo` otherwise bash will think 'originally' is a command:
`cat filt_vs_original | sed 's,^orgin,echo "origin,g' | sed 's,filtering_results,filtering_results",g' > compare_filtered`    

Then launch the job (counting the number of lines in each file takes a long time, otherwise we could just do that on the head node):

`launcher_creator.py -n compare_filtered -a Sailfin_RNASeq -q normal -t 0:10:00 -e lukereding@utexas.edu -j compare_filtered`    
`qsub launcher.sge`

Actually, it saved in one of the job files, so run: `grep ^originally compare_filtered.o2953770 > compare_filtered`     

>originally: 35832124 , filtered: 32671630 >> filtering_results         
originally: 42569481 , filtered: 38775246 >> filtering_results           
originally: 49402160 , filtered: 45174284 >> filtering_results         
originally: 47989509 , filtered: 43933380 >> filtering_results            
originally: 58941376 , filtered: 53920567 >> filtering_results            
originally: 66531156 , filtered: 60649980 >> filtering_results          
originally: 59410809 , filtered: 54536624 >> filtering_results         
originally: 56965320 , filtered: 52016187 >> filtering_results         
originally: 64199757 , filtered: 59026770 >> filtering_results        
originally: 66090458 , filtered: 60478621 >> filtering_results          
originally: 68997077 , filtered: 63162381 >> filtering_results          


### adaptors
Adaptors at at the end of the read and anchor the read to the flow cell. There is likely little adpator contamination, so I skipped this step.

---------------

## getting a transcriptome

I first decided to use the guppy transcriptome:       
`wget ftp://ftp.tuebingen.mpg.de/ebio/publication_data/esharma/guppy_trans/trin_cuff_v14_cdhit90.fa.gz`   # download the transcriptome     
`gunzip trin_cuff_v14_cdhit90.fa.gz` #unzip it     
`mv trin_cuff_v14_cdhit90.fa guppy_transcriptome` # rename it    


## making transcriptome annotations
### doing the blasts
 The transcriptome is just a list of sequences that are likely genes; we don't know the identites of the genes. For this, we have to do a series of BLASTs, blasting each sequence against know gene sequences and names our genes based on the results. The Swiss-Prot database contains a bunch of sequences that we can blast against (I think). Note that I'm transitioning here from using the walkthrough provided [here](https://wikis.utexas.edu/display/bioiteam/Introduction+to+RNA+Seq+Course+2015) to Misha's workflow, which can be found [here](https://github.com/z0on/annotatingTranscriptomes/blob/6751534ed87b0893242be97adfd00de1012a08a3/annotating%20trascriptome.txt)


`echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" > get_prot`        
`launcher_creator.py -j get_prot -l get_prot -a Sailfin_RNASeq -q normal -t 1:00:00 -e lukereding@utexas.edu`          
`qsub get_prot`           

Get the annotations next:     
`echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz" > get_go`       
`launcher_creator.py -j get_go -l get_go -a Sailfin_RNASeq -q normal -t 4:00:00 -e lukereding@utexas.edu`
`qsub get_go`


For some reason this didn't work for me. TACC threw a bunch of errors into the output of STDERR. Instead, I downloaded each to my local machine. The first file, `uniprot_sprot.fasta.gz`, downloaded in like 10 seconds. The second file took about three minutes. Unzip each file with `gunzip`. `uniprot_sprot.fasta` is ~250 MB (a few seconds to transfer to TACC via `scp`) and the annotation file is almost 11 GB (~5 min to transfer to TACC via `scp`).        

Next we have to index the uniprot fasta file:

`module load blast` # load the blast module
`echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" > make_index` # create job script         
`launcher_creator.py -j make_index -n make_index -l index -a Sailfin_RNASeq -e lukereding@utexas.edu` # create launcher file          
`qsub index` # submit job           

The blast command looks something like this:

`blastx -query guppy_transcriptome -db uniprot_sprot.fasta -evalue 0.0001 -num_threads 36 -num_descriptions 5 -num_alignments 5 -out guppy_transcriptome.br`        

What do the arguments mean?      
- query: the transcriptome you want to annotate. Misha does this by breaking up the transcriptome into a bunch of different pieces. I'm going to first try without doing this.       
- db: the database to blast against, i.e., the swiss prot database       
- evalue: e-value threshold for saving hits. Above is the same e-value that both Misha uses and the Fraser et al. paper uses.       
- num_threads: number of threads or cores to use. default is 1       
- num_descriptions: Number of database sequences to show one-line descriptions for (??)        
- out: name of output file        

I am going to try submitting this job for 24 hours to three nodes (note that each node has 12 cores, 12 * 3 = 36, the number of threads in our call to `blastx` above.)

`echo "blastx -query guppy_transcriptome -db uniprot_sprot.fasta -evalue 0.0001 -num_threads 36 -num_descriptions 5 -num_alignments 5 -out guppy_transcriptome.br" > blast` # make job script              
`launcher_creator.py -j blast -n blast_job -a Sailfin_RNASeq -e lukereding@utexas.edu -t 23:55:00 -q normal` # create launcher file. note that this must do to the "normal" quene as the development queue has a time limit of 1 hour                
`cat launcher.sge | perl -pe 's/12way .+$/1way 36/' >blast_job_36`      # substitute the wayness           
`qsub blast_job_36` # send it in           

This took roughly 8 hours to run.     
----------------

### an aside

I find it's always good practice to check the number of things (e.g. genes) in the files you're working with for datasets this large. Note that `grep "^>" guppy_transcriptome | wc -l` gives `74567`, which seems like far too many genes. Of course, this number might be so high because we're using a transcriptome and therefore splice varients will be represented by different transcripts. This number still seems too high to me though; I'll need to go back to the original paper and make sure it's correct.      

Below, we make some files we need to be able to use one of Misha's many scripts. Doing that involves naming all the genes in our transcriptome. Oddly, when I go through this process and call `wc -l transcriptome_seq2iso.tab` to give the number of genes in this file, the result is `18100`, which seems much more reasonable to me. I haven't been able to figure out why there is such a large discrepancy here.

--------------------
### processing blast results
According to Misha's, script: "if there are no components or isogroups in your transcriptome, create fake isogroup designations (=contigs)". The alternative to this is to use Trinity to assemble a transcriptome _de novo_, which we could do but don't want to do. I ran these lines according to his script:      

>Don't run these lines!          
`grep ">" guppy_transcriptome.fasta | perl -pe 's/>(\S+)\s.+/$1\tisogroup$1/' >transcriptome_seq2iso.tab`            
`cat guppy_transcriptome.fasta | perl -pe 's/>(\S+).+/>$1 gene=isogroup$1/' >transcriptome_iso.fasta`             


These lines don't work. The goal: create a `transcriptome_seq2iso.tab` file that has one column representing trasncript ID and a second column with gene ID. I'm going to assume that 1 transcript = 1 gene. Then we want to create a `transcriptome_iso.fasta` file in which the headers are the same transcript ID, tab, then 'gene=geneID'. These must be the same as the names in the `transcriptome_seq2iso.tab` file. Finally, these names must match the 'Query = ' lines in the blast search results (the `.br` file that was generated from blast).

`grep ">" guppy_transcriptome | cut -d"_" -f2,3 | sed 's,_m,,g' | sed 's,\([0-9]*\)\.\([0-9]*\),comp\1\2\tisogroup\1\2,' > transcriptome_seq2iso.tab`          


Next we need to change the headers on the transcriptome. The new header is tab delimited two column of transcript ID and gene ID in the format 'gene=geneID'.       

`cat guppy_transcriptome | sed 's,^>CUFF_\([0-9]*\),>comp\1,g'  | sed 's,_m,,g' | sed 's,comp\([0-9]*\).\([0-9]*\),comp\1\2 gene=isogroup\1\2,g' > transcriptome_iso.fasta`          

Finally, we change the blast search results; I save them under a new name, `guppy_blast_results_corrected.br`:          

`cat *.br | sed 's,^Query= CUFF_\([0-9]*\),Query= \1,g' | sed 's,_m.,,g' | sed 's,^Query= \([
0-9]*\),Query= comp\1 gene=isogroup\1,g' > guppy_blast_results_corrected.br`          


cat *.br | sed 's,^Query= CUFF_\([0-9]*\),Query= \1,g' | sed 's,_m.,,g' | sed 's,^Query= \([0-9]*\),Query=  gene=isogroup\1,g' > test.br


What do these files looks like now?           


`head transcriptome_seq2iso.tab`                 
>comp10004275477	isogroup10004275477       
comp10006275362	isogroup10006275362
comp10006275360	isogroup10006275360
comp10009275364	isogroup10009275364
comp10013275481	isogroup10013275481
comp10015275438	isogroup10015275438
comp10016275544	isogroup10016275544
comp10017275428	isogroup10017275428
comp10018327682	isogroup10018327682
comp10017252	isogroup10017252

`head transcriptome_iso.fasta`             
>comp10004275477 gene=isogroup10004275477
ATGTGCACACCTGATGTTGCTCACAGGGGTCCTCGGTGTGCTCAGGGAGCTTGTGTATAC
GCCCAGACGCATGAAAAATTAGAGGGAACATTGGTCGGTTCGCTTCTCTTTTCCTCCACA
CCTGAACTTACCTGTCCGCCCTCCACCTCAGCCTCCACATCACCTACAAAGAACTCACTC
CATCAAGACCTCACACTGGAACCAGGACCACTTCAAGGAGACCCACCAAGGGATCATACT
TTGGGAACCACTGTATTACCCCCCCGGCGGAAGTCTATCTTCATATCAAGTAAAATCAAT
TTAATTTTTCCTGAAGTTCCATTTTCCATAAAGAGTAAAATTATTGTGATGGAATGA             
>comp10006275362 gene=isogroup10006275362
ATGAAGGATGCAAGCAACAAAAGGGTGAGCAGATTGGTTAACGTGAAGCCCCTCACACCC
TTTATGGTGCCAGATGTCATGTTGGCTTGGTGTCCAGCTCAGAAACAGGATCAGGTAACT

`head -25 guppy_blast_results_corrected.br | tail -10`           
>Query= comp10004275477 gene=isogroup10004275477
>
>Length=357
>
>
>***** No hits found *****
>
>
>
>Lambda      K        H        a         alpha



Because we are using a previously assembled transcriptome, we don't need to use _cd hit_ to estimate isogroups (genes), the next step in his pipeline.     

Next, we want to use the results from our blast to extract the gene names for each of our transcripts. `getGeneNameFromUniProtKB.pl` is a Perl script that Misha wrote and is available [here](https://github.com/z0on/annotatingTranscriptomes/blob/master/getGeneNameFromUniProtKB.pl). The easiest way is to download the repository from github and use `scp` to transfer it to TACC. Note that you may need to do `chmod +x getGeneNameFromUniProtKB.pl` to change permissions to ensure you can execute the script. Note also that you need to include the word `perl` before the name of the script; `which perl` yeilds `/opt/apps/perl/5.14/bin/perl`, not `/usr/bin/perl/` which in the shebang of Misha's script.

`perl getGeneNameFromUniProtKB.pl blast=guppy_blast_results_corrected.br prefix=guppy fastaQuery=transcriptome_iso.fasta`

`head guppy_iso2hit.tab`
>sp|Q92817|EVPL_HUMAN	1e-180
isogroup20402102620	sp|P26990|ARF6_CHICK	3e-127
isogroup50792238338	sp|P35072|TCB1_CAEBR	3e-15
isogroup18699303674	sp|A1L390|PKHG3_HUMAN	3e-178
isogroup3307930504	sp|P11279|LAMP1_HUMAN	3e-61
isogroup1679495026	sp|P34205|PHR_CARAU	1e-13
isogroup1681895243	sp|Q92536|YLAT2_HUMAN	7e-76
isogroup36729173669	sp|Q91453|STXB_SYNHO	6e-22
isogroup18931304870	sp|Q8QZR1|ATTY_MOUSE	1e-20
isogroup2834526737	sp|Q04799|FMO5_RABIT	2e-100


`head guppy_iso2gene.tab`
>	Envoplakin OS=Homo sapiens GN=EVPL PE=1 SV=3 E(blastx)=1e-180
isogroup20402102620	ADP-ribosylation factor 6 OS=Gallus gallus GN=ARF6 PE=2 SV=3 E(blastx)=3e-127
isogroup50792238338	Transposable element Tcb1 transposase OS=Caenorhabditis briggsae PE=3 SV=1 E(blastx)=3e-15
isogroup18699303674	Pleckstrin homology domain-containing family G member 3 OS=Homo sapiens GN=PLEKHG3 PE=1 SV=1 E(blastx)=3e-178
isogroup3307930504	Lysosome-associated membrane glycoprotein 1 OS=Homo sapiens GN=LAMP1 PE=1 SV=3 E(blastx)=3e-61
isogroup1679495026	Deoxyribodipyrimidine photo-lyase OS=Carassius auratus GN=phr PE=2 SV=1 E(blastx)=1e-13
isogroup1681895243	Y+L amino acid transporter 2 OS=Homo sapiens GN=SLC7A6 PE=1 SV=3 E(blastx)=7e-76
isogroup36729173669	Stonustoxin subunit beta OS=Synanceia horrida PE=1 SV=3 E(blastx)=6e-22
isogroup18931304870	Tyrosine aminotransferase OS=Mus musculus GN=Tat PE=1 SV=1 E(blastx)=1e-20
isogroup2834526737	Dimethylaniline monooxygenase [N-oxide-forming] 5 OS=Oryctolagus cuniculus GN=FMO5 PE=1 SV=2 E(blastx)=2e-100

Next we want to extract to GO terms associated with each gene:

`perl getGOfromUniProtKB.pl blast=guppy_blast_results_corrected.br prefix=transcriptome fastaQuery=transcriptome_iso.fasta`

Let's look at the file that results:          

`head transcriptome_iso2go.tab`

>	GO:0001533;GO:0005737;GO:0030057;GO:0070062;GO:0045111;GO:0030674;GO:0005198;GO:0008544;GO:0031424;GO:0030216;GO:0018149
isogroup20402102620	GO:0005938;GO:0005737;GO:0030139;GO:0005768;GO:0070062;GO:0031527;GO:0005925;GO:0005794;GO:0043209;GO:0005886;GO:0055038;GO:0001726;GO:0005525;GO:0030866;GO:0090162;GO:0097284;GO:0001889;GO:0033028;GO:0030838;GO:0090004;GO:0034394;GO:0036010;GO:0015031;GO:0060998;GO:0051489;GO:0035020;GO:0031529;GO:0007264;GO:0016192
isogroup50792238338	GO:0005634;GO:0003677;GO:0015074;GO:0006313
isogroup18699303674	GO:0005089;GO:0035023
isogroup3307930504	GO:0097208;GO:0044194;GO:0005737;GO:0030425;GO:0010008;GO:0009897;GO:0070062;GO:0005887;GO:0005770;GO:0005764;GO:0042470;GO:0016020;GO:0005771;GO:0043025;GO:0048471;GO:0061474;GO:0005886;GO:0042383;GO:0008021;GO:0019899;GO:0001618;GO:0048102;GO:0006914;GO:0072594;GO:0090160;GO:0008626;GO:0043323;GO:0045954;GO:0050821;GO:1902513
isogroup1679495026	GO:0003904;GO:0003677;GO:0006281;GO:0018298
isogroup1681895243	GO:0016323;GO:0005887;GO:0005886;GO:0015171;GO:0015297;GO:0015179;GO:0003333;GO:0006865;GO:0007596;GO:0006520;GO:0006811;GO:1902475;GO:0015807;GO:0050900;GO:0006461;GO:0055085;GO:0006810
isogroup36729173669	GO:0005576;GO:0044179
isogroup18931304870	GO:0005739;GO:0016597;GO:0080130;GO:0004838;GO:0030170;GO:0006103;GO:0009058;GO:0006536;GO:0006559;GO:0051384;GO:0046689;GO:0006979;GO:0006572
isogroup2834526737	GO:0005789;GO:0016021;GO:0050660;GO:0004499;GO:0050661


------------------
## mapping pt. 1: BWA

Meanwhile, I should have been trying to map the reads to the guppy genome. I'll follow the example on the RNA seq walkthrough and use BWA to map. Let's first index the guppy transcriptome:

`module load bwa/0.7.7`        
`echo "bwa index -a bwtsw transcriptome_iso.fasta" > index_guppy`
`launcher_creator.py -j index_guppy -n index_guppy_job -a Sailfin_RNASeq -e lukereding@utexas.edu -t 1:00:00`           

Now let's make the mapping commands using BWA-MEM:    

Create job file.       
`for file in *.filtered; do echo "bwa mem transcriptome_iso.fasta $file > ${file%.fastq.filtered}.sam" >> mapping_commands; done`           

Create launcher.        
`launcher_creator.py -j mapping_commands -n mapping_commands_job -l mapping_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 23:50:00`

Substitute wayness so that things run faster:          
`cat mapping_commands_job | perl -pe 's/12way .+$/4way 48/' >mapping_commands_jobscript`              
Submit the job:         
`qsub mapping_commands_jobscript`            

---------------
### SAM conversion, sorting, and indexing

We now convert the SAM files that resulted from the previous step to BAM files, which are less unwieldy.

`for file in SRR*.sam; do echo "samtools view -b -S $file > ${file%.sam}.bam && samtools sort ${file%.sam}.bam ${file%.sam}.bam.sort && samtools index ${file%.sam}.bam" >> bam_commands; done`

Launch it:     
`launcher_creator.py -j bam_commands -n bam_commands_job -l bam_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 2:00:00`
`qsub bam_commands_job`


### assessing mapping quality

Now we can assess mapping quality for the reads from each individual like:

`samtools flagstat SRR1165201_1.bam`

To get an overall sense of the percentage of reads that mapped, I ran a job with the command `samtools flagstat SRR*.bam`. The result (at the bottom of `BAM_job.o2959735`), shows `24006209 + 0 mapped (50.04%:nan%)`, suggesting that only half the reads mapped. Not great.

-----------------

## gene counting


`module load samtools`          

For each individual (i.e. each BAM file), we can use `idxstats` in `samtools` to extract a list of genes and their counts:
`samtools view SRR1166372_1.bam | cut -f3 | sort| uniq -c | sort -k1nr | cut -f1-2`

Let's do this for each individual:

`for file in *.bam.sort.bam; do echo "samtools view $file | cut -f3 | sort| uniq -c | sort -k1nr | cut -f1-2 > ${file%.bam.sort.bam}.counts" >> counts_job; done`            
`launcher_creator.py -j counts_job -n counts_job_job -l counts_job_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 5:00:00`              
`qsub counts_job_job`

**note:** also see the bash script [here](https://www.biostars.org/p/14531/) for a possible, better alternative. also see notes [here](http://davetang.org/wiki/tiki-index.php?page=SAMTools#Basic_usage)   

-----------------

Because of the terrible mapping success we had, I'm going to try a different approach and use Bowtie2 instead of BWA to map the reads. This sort of compromises the nice linear structure I've kept up until this point, but I think it's necessary to document the steps I've taken.

-----------------

## mapping pt. 2: Bowtie2

Now I retry mapping, this time with Bowtie2 to see if more of our reads will map to the guppy transcriptome.

To keep this attempt separate from the BWA attempt, I'll contain all the Bowtie mapping files to a folder within `molly_arts` on TACC.

`mkdir bowtie2_mapping && cd bowtie2_mapping`
Copy all filtered fastq files:
`echo "find .. -name '*.filtered' -exec cp {} . \;" > copy`
`launcher_creator.py -j copy -n copy_job -l copy_job -a Sailfin_RNASeq -e lukereding@utexas.edu -t 1:00:00`
Rename so that the filtered files end in `.fasta` otherwise bowtie will throw an error:
`for file in *.filtered; do mv $file $file.fasta; done`

Before mapping, create an index of the `guppy_transcriptome` FASTA file:

`module load samtools`
`module load bowtie/2.1.0`
`bowtie2-build ../guppy_transcriptome guppy_transcriptome`


Before submitting the job, let's look at the main function call and try to understand the arguments:

`bowtie2 --local --very-fast-local -f -k 5 -x ../guppy_transcriptome.fai -U ../XX.filtered -S XX.sam`

* --local : From the [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml):

>When the --local option is specified, Bowtie 2 performs local read alignment. In this mode, Bowtie 2 might "trim" or "clip" some read characters from one or both ends of the alignment if doing so maximizes the alignment score.       

* --very-fast-local : not sure
* -f : tells bowtie2 that these are FASTA files
* -k 1 : again, from the manual:

>In -k mode, Bowtie 2 searches for up to N distinct, valid alignments for each read, where N equals the integer specified with the -k parameter. That is, if -k 2 is specified, Bowtie 2 will search for at most 2 distinct alignments. It reports all alignments found, in descending order by alignment score.

>Note that Misha uses `-k 5`

* -x : the index of the reference file
* -U : reads to be aligned
* -S : file to write SAM alignments to


Let's write the job file:
`module load bowtie/2.1.0`
`for file in *.filtered.fasta; do echo "bowtie2 --local --very-fast-local -k 5 -x guppy_transcriptome -U $file -S ./${file%.fastq.filtered.fasta}.sam" >> bowtie2_mapping_commands; done`

Create launcher.        
`launcher_creator.py -j bowtie2_mapping_commands -n bowtie2_mapping_commands_job -l bowtie2_mapping_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 23:50:00`

Substitute wayness so that things run faster:          
`cat bowtie2_mapping_commands_job | perl -pe 's/12way .+$/4way 48/' >mapping_commands_jobscript`              
Submit the job:         
`qsub mapping_commands_jobscript`            

This took one hour and fifteeen minutes to run.

I now follow the same steps as above, repeated here for completeness.

#### SAM conversion, sorting, and indexing

We now convert the SAM files that resulted from the previous step to BAM files, which are less unwieldy.

`for file in *.sam; do echo "samtools view -b -S $file > ${file%.sam}.bam && samtools sort ${file%.sam}.bam ${file%.sam}.bam.sort && samtools index ${file%.sam}.bam" >> bam_commands; done`

Launch it:     
`launcher_creator.py -j bam_commands -n bam_commands_job -l bam_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 2:00:00`
`qsub bam_commands_job`

For whatever reason, the indexing of the SAM files didn't work. I'll do a separate job to try that:
`cat bam_commands | sed 's,^\(.*\) && \(samtools index .*\),\2,g' > bam_index`
This puts all the indexing commands in one file. Then we have have index the sorted BAM files, not the unsorted ones (I wrote the initial line wrong):
`cat bam_index | sed 's,\(S.*\).bam,\1.bam.sort.bam,g' > bam_index_sort`
`launcher_creator.py -j bam_index_sort -n bam_index_sort_job -l bam_index_sort_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 2:00:00`
`qsub bam_index_sort_job`



#### assessing mapping quality

Now we can assess mapping quality for the reads from each individual like:

`samtools flagstat SRR1165201_1.bam`

To get an overall sense of the percentage of reads that mapped, I ran a job with the command `samtools flagstat SRR*.bam`. The result (at the bottom of `BAM_job.o2959735`), shows `25424751 + 0 mapped (49.15%:nan%)`, suggesting that, again, only half the reads mapped. Not great, but at least we're getting some sort of consistent result between the two mapping algorithms.


---------------------

## analysis of gene counts from BWA

At this point, I have series of files, one per individual, containing the number of mRNA reads mapped to get gene. The goal is to generate some sort of spreadsheet where the rows represent all possible genes and the column represent individuals. The cells are filled with the number of transcripts mapped to a given gene for a given individual.

I wrote a script, `merge_counts.R`, to automate this. The script should be in this directory / repository. The script is run from the command line like `Rscript merge_counts.R /path/to/directory`. The script will look in `path/to/directory` for all files that end in `.counts`. These files come from the above *gene counting* step. These are files with two columns separated with whitespace (so that `read.table` in R can read them). The first column are counts, the second column is gene names. The script merges all these `counts` together, naming the columns of the resulting dataframe based on the names of the `.counts` files. The resulting dataframe is saved as `counts.csv` in `/path/to/directory`. It has n+1 columns, where n = number of `.counts` file (or individuals). If no gene is expressed in one inidividual but is expressed in other individuals, I change the NA supplied by the `merge` function in R to a zero, meaning that no expression of that gene was detected. You can run this R script on TACC and then move the resulting csv file to your local machine.





### another aside: useful samtools commands

*Taken from [this](http://davetang.org/wiki/tiki-index.php?page=SAMTools#Basic_usage) website.*

Covert a SAM file to a sorted BAM file in one step:

`samtools view -bS file.sam | samtools sort - file_sorted`

Filter out unmapped reads from a BAM file, save to a new BAM file:

`samtools view -h -F 4 -b blah.bam > blah_only_mapped.bam`

Create a FASTQ file from a BAM file (only aligned reads):

`bam2fastq -o blah_unaligned.fastq --no-unaligned blah.bam`

Fastest way to count number of reads in a BAM file:

`samtools idxstats in.bam | awk '{s+=$3+$4} END {print s}`   # number of reads            
`samtools idxstats in.bam | awk '{s+=$3} END {print s}'    `#number of mapped reads

--------------

# Redoing everything with the Amazon molly genome.

I am convinced that the guppy transcriptome is a poor reference to map to. Not only am I getting only ~50% of my reads mapped, regardless of whether I use Bowtie or BWA, but a previous paper by the same group that initially characterized the guppy transcriptome found that 40% of their reads from sailfins mapped to the guppy transcriptome, which is largely consistent with my findings, suggesting I'm not doing something terribly wrong.

I am going to re-do everything use the Amazon molly genome as a reference.

Make a new folder called `amazon` to keep things tidy:

`mkdir amazon`        
`cd amazon`      

Copy the original fastq files int the new folder so that everything is self-contained:

`cp ../*.fastq.filtered .`

Download the amazon molly reference genome from emsembl:

`curl -O ftp://ftp.ensembl.org/pub/release-83/fasta/poecilia_formosa/dna/Poecilia_formosa.PoeFor_5.1.2.dna.toplevel.fa.gz`        
`gunzip Poecilia_formosa.PoeFor_5.1.2.dna.toplevel.fa.gz`    

This results in `Poecilia_formosa.PoeFor_5.1.2.dna.toplevel.fa`. I'm going to rename it to make it easier to work with:

`mv Poecilia_formosa.PoeFor_5.1.2.dna.toplevel.fa amazon_genome.fasta`

Let's also pull the GTF file from ensembl while we're at it:

`curl -O ftp://ftp.ensembl.org/pub/release-83/gtf/poecilia_formosa/Poecilia_formosa.PoeFor_5.1.2.83.gtf.gz`     
`gunzip *.gz`    


---------------

### mapping with BWA and the amazon molly genome

We first need to index the amazon molly genome prior to mapping:

`module load bwa/0.7.7`        
`echo "bwa index -a bwtsw amazon_genome.fasta" > index_amazon`
`launcher_creator.py -j index_amazon -n index_amazon_job -l index_amazon_job -a Sailfin_RNASeq -e lukereding@utexas.edu -t 1:00:00`           





Now let's make the mapping commands using BWA-MEM:    

Create job file.       
`for file in *.filtered; do echo "bwa mem amazon_genome.fasta $file > ${file%.fastq.filtered}.sam" >> mapping_commands; done`           

Create launcher.        
`launcher_creator.py -j mapping_commands -n mapping_commands_job -l mapping_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 23:50:00`

Substitute wayness so that things run faster:          
`cat mapping_commands_job | perl -pe 's/12way .+$/4way 48/' >mapping_commands_jobscript`              
Submit the job:         
`qsub mapping_commands_jobscript`           


#### SAM conversion, sorting, and indexing

We now convert the SAM files that resulted from the previous step to BAM files, which are less unwieldy.

`for file in SRR*.sam; do echo "samtools view -b -S $file > ${file%.sam}.bam && samtools sort ${file%.sam}.bam ${file%.sam}.bam.sort && samtools index ${file%.sam}.bam" >> bam_commands; done`      

`launcher_creator.py -j bam_commands -n bam_commands_job -l bam_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 2:00:00`   

`qsub bam_commands_job`    

Oops:

`cat bam_commands | sed 's,\(.*\)\(samtools index .*\),\2,g' > bam_indexing`    

`launcher_creator.py -j bam_indexing -n bam_indexing_job -l bam_indexing_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 1:00:00`  

`qsub bam_indexing_job`

Oops again:

`cat bam_indexing | sed 's,.bam,.bam.sort.bam,g' > bam_index`

`launcher_creator.py -j bam_index -n bam_index_job -l bam_index_job -a Sailfin_RNASeq -e lukereding@utexas.edu -t 0:15:00`

`qsub bam_index_job`

#### assessing mapping quality

`echo "samtools flagstat SRR*.bam > mapping_results" > mapping_quality`     

`launcher_creator.py -j mapping_quality -n mapping_quality_job -l mapping_quality_job -a Sailfin_RNASeq -e lukereding@utexas.edu`   

`qsub mapping_quality_job`

`cat`ing `mapping_results` yeilds `50878768 + 0 mapped (95.40%:nan%)`: 95% of the reads mapped!

### counting genes

#### using samtools to count genes

Instead of doing what I did above, a (likely) better way of doing things is to use `samtools idxstats bam_file`. This returns a list with four tab-delimited columns. The first is the gene/transcript, the sequence length, the number of reads mapped, and the number of reads unmapped.

`for file in *.bam.sort.bam; do echo "samtools idxstats $file > ${file%.bam.sort.bam}.counts" >> count_genes; done`

`launcher_creator.py -j count_genes -n count_genes_job -l count_genes_job -t 0:15:00 -a Sailfin_RNASeq -e lukereding@utexas.edu`

`qsub count_genes_job`

In reality, this isn't super useful. The `.counts` files look like this:

>KI519610.1	7371162	574071	0         
KI519611.1	6947489	609704	0            
KI519612.1	6918293	628381	0          
KI519613.1	6837434	721852	0       
KI519614.1	6763396	389486	0       
KI519615.1	6332452	361961	0         
KI519616.1	6276451	577374	0          

The first column is the chromosome and the next two of the regions between which a given feature occurs. The last column is the number of reads mapping to that feature. We want gene names, not locations on chromosomes, so counting things this way isn't super useful.

~~Before I downloaded the .gtf file; instead, to follow along with the intructions on TACC, I'm going to download the .gff file:~~

~~`curl -O ftp://ftp.ensembl.org/pub/release-83/gff3/poecilia_formosa/Poecilia_formosa.PoeFor_5.1.2.83.gff3.gz`~~

~~`gunzip Poecilia_formosa.PoeFor_5.1.2.83.gff3.gz`~~

~~Rename:~~

~~`mv Poecilia_formosa.PoeFor_5.1.2.83.gff3 amazon.gff3`~~

#### using bedtools to count genes

I'll use bedtools to count genes:

`module load bedtools`

`for file in *.bam.sort.bam; do echo "bedtools multicov -bams $file -bed *.gtf > ${file%.bam.sort.bam}.gff" >> counting; done`

`launcher_creator.py -j counting -n counting_job -l counting_job -a Sailfin_RNASeq -e lukereding@utexas.edu`

 `qsub counting_job`


This results in a bunch of `.gff` files that look like this:

> KI519610.1	ensembl	gene	772	45558	.	+	.	gene_id "ENSPFOG00000016007"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding";	3098
KI519610.1	ensembl	transcript	772	45558	.	+	.	gene_id "ENSPFOG00000016007"; gene_version "1"; transcript_id "ENSPFOT00000016137"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";	3098
tacc:/scratch/02535/reding/molly_arts/amazon$ head -2 SRR1166368_1.gff | grep "gene_source"
KI519610.1	ensembl	gene	772	45558	.	+	.	gene_id "ENSPFOG00000016007"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding";	3098
KI519610.1	ensembl	transcript	772	45558	.	+	.	gene_id "ENSPFOG00000016007"; gene_version "1"; transcript_id "ENSPFOT00000016137"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";	3098

The last number on each line is the number of reads that map to that feature.

The problem here is that I want to just focus on genes, and the `.gff` files contains exons, UTRs, etc. that I don't care about.

I thought that maybe I could just `grep` all line that contain `gene_name`, but `grep "^.*" SRR1166368_1.gff | wc -l` yeilds __831661__ (i.e. there are 831,661 lines in the gff file) while `grep "gene_name" SRR1166368_1.gff | wc -l` yeilds __742961__ (way too many to represent different genes). Comparing `grep "gene_name" SRR1166368_1.gff | wc -l` to the lines that __don't__ contain `gene_name` (using inverse grep, i.e., `grep -v "gene_name" SRR1166368_1.gff | wc -l`) doesn't show patterns: exons, UTRs, CDS, etc. are represented in both. I don't know why some lines contain the gene name and other don't.

Note that if we use `awk`--`cat P*.gtf | awk '$3=="gene"' | wc -l`--we get __24,354__, which is probably pretty close to the number of genes. We can use the same `awk` command on our `.gff` files to extract only the lines containing genes and not more specific stuff that we might not want. (For a primer on how to read GTF files, look [here](http://useast.ensembl.org/info/website/upload/gff.html?redirect=no))

Before we do that though, I'm going to re-do the call to `bedtools multicov` above. From the [bedtools docs] (http://bedtools.readthedocs.org/en/latest/content/tools/multicov.html), if you feed `multicov` more than one sorted bam file, it'll just add columns onto the resulting `.gff`, one column per `.bam` file. This will make it easier to wrangle the data for use is DESeq2 as long as we are careful keep track of what sample goes with what column.


`echo "bedtools multicov -bams SRR1161450_1.bam.sort.bam SRR1165203_1.bam.sort.bam SRR1166368_1.bam.sort.bam SRR1166371_1.bam.sort.bam SRR1161451_1.bam.sort.bam SRR1166366_1.bam.sort.bam SRR1166369_1.bam.sort.bam SRR1166372_1.bam.sort.bam SRR1165201_1.bam.sort.bam SRR1166367_1.bam.sort.bam SRR1166370_1.bam.sort.bam -bed P*.gtf > results.gff" > bedtools_commands`

`launcher_creator.py -j bedtools_commands -n bedtools_commands_job -l bedtools_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 10:00:00`

`qsub bedtools_commands_job`

__note__: this took 6.5 hours to run for some reason.

From this we have `results.gff`. To make our lives simple, we wish to only extract features that represent genes. We also want to extract the ensembl gene IDs, gene names, and counts. Note that for some reason the number of columns on each line can be difference depending on whether 'gene_name' or 'gene_source' is in the third column (using ; as a delimiter):


`cat results.gff | awk '$3=="gene"' | cut -f9- | cut -d';' -f1,3,4,5- | cut -d' ' -f2,4- | cut -d ';' -f1,2,4- | awk -F ";" '{print $1,$2,$(NF)}' | sed 's,\",,g' > counts.tab`


> __Explanation:__
`cat` is pushing the `result.gff` file to stdout        
`awk` is selecting lines in which the third column (tab-delimited) says 'gene'          
`cut -f9-` is getting rid of the first 8 columns, which contain the chromosome and other stuff we don't really care about         
`cut -d';' -f1,3,4,5-` is changing the delimiter for `cut` from tabs (the default) to semi-colons. It then selects the 1st, 3rd, and 5th - end columns, which contain the gene ID, gene name, and counts for each individual           
`cut -d' ' -f2,4-` changes the delimiter to a space and selects only the gene id itself (getting rid of the 'gene_id' label), the gene name itself, and the counts           
`awk -F ";" '{print $1,$2,$(NF)}'` is changing the default delimiter for `awk` from a tab to a semicolon, then printing the first, second, and last columns. Some of the rows contain an extra variable in the third column, and this part of the one-liner gets rid of it             
`sed 's,\",,g'` gets rid of all double quotation marks           


So we've taken a file with lines that look like this: 

> KI519610.1	ensembl	gene	772	45558	.	+	.	gene_id "ENSPFOG00000016007"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding";	2529	2920	3098	2600	1778	2569	1948	962	1849	2398	1511

and turned it into something much more useful to us:

> ENSPFOG00000016007  ensembl 	2529	2920	3098	2600	1778	2569	1948	962	1849	2398	1511

`wc -l counts.tab` yeilds `24354`. 

We now can transfer this file to our local machine, using `scp` from our local machine like

`scp reding@lonestar.tacc.utexas.edu:/scratch/02535/reding/molly_arts/amazon/counts.tab /Users/lukereding/Desktop/molly_arts/`

-------------

### creating annotation files to use for GO analyses


The steps here are repeats of those listed above for getting GO terms from the guppy transcriptome. They are repeated here in the sake of completeness.

Get the annotations:
`wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`         
`wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`

`gunzip *.gz`

~~`module load blast`~~      
~~`echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" > make_index`
`launcher_creator.py -j make_index -n make_index -l index -a Sailfin_RNASeq -e lukereding@utexas.edu`
`qsub index`~~


~~`echo "blastx -query amazon_genome.fasta -db uniprot_sprot.fasta -evalue 0.0001 -num_threads 36 -num_descriptions 5 -num_alignments 5 -out amazon_genome.br" > blast`                     
`launcher_creator.py -j blast -n blast_job -a Sailfin_RNASeq -e lukereding@utexas.edu -t 23:55:00 -q normal`               
`cat launcher.sge | perl -pe 's/12way .+$/1way 36/' >blast_job_36`          
`qsub blast_job_36~~


This keeps timing out on me. We will use Misha's method of breaking up the genome into little chunks, running the chunks separately, then putting concatenating the results together.

I'll use Misha's script `splitFasta.pl` to split it into 

`perl splitFasta.pl amazon_genome.fasta 40`     
`module load blast`
`ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 3 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl`      
`launcher_creator.py -j bl -n split_blast_job -t 1:00:00 -a Sailfin_RNASeq -e lukereding@utexas.edu`        
# change wayness
`cat launcher.sge | sed 's,-pe 12way 48,-pe 4way 120,' > blast_job'`
`qsub blast_job`

-----------------


At this point, I'm going to transition from `bash` to `R` and record what I'm doing in the `deseq_sailfin_analysis.Rmd` file in this repo.



---------------

### from the supplemental materials:
short read archive -- SRP036185

_*for reference:*_         
_(a) RNA-seq data collection and analysis RNA was extracted from whole brain tissue using RNeasy Lipid Tissue Mini Kit (Qiagen Hilden, Germany). We prepared separate RNA sequencing libraries from whole brains for each individual, using unique index sequences from the Illumina Tru-Seq RNA kit following manufacturers instructions. Sequencing libraries were constructed and sequenced on three lanes of an Illumina HiSeq 2000 at the HudsonAlpha Genomic Services Laboratory (Huntsville, Alabama) in April 2012. The reference P. reticulata assembly was constructed from a data set containing > 450 million 100-bp paired end reads, which were filtered for high quality sequence and normalized in-silico to compress the range in kmer abundance. We used SeqMan NGEN 4.1.2 (Madison, WI) [74,75] to perform the assembly. Contigs from the assembly were annotated by blastx queries against SwissProt (database downloaded Oct 6, 2012), UniProt/Trembl (Nov 28, 2012), and nr (Dec 11, 2012). Default parameters were used in the blastx queries, with e-value cutoff of 10-4 . The assembly and individual reads will be deposited in a publicly-available archive before publication. Reads were mapped to the reference assembly using Bowtie 2 v 2.0.0 on a server running Red Hat Enterprise Linux 6.1. We used a seed size of 20 bp, with no mismatches allowed in the seed (run options: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75). We retained mappings with quality scores > 30 (< 0.001 probability that the read maps elsewhere in the reference), and applied an abundance filter to retain only transcripts represented by more than 1 count per million reads in at least three samples. This filtering resulted in 31,869 transcripts remaining in the data set. We used the number of reads mapping to each of those transcripts, along with TMM-normalized library sizes [76] to analyse differential expression. We obtained 462,537,724 100-bp reads that passed the machine quality filter, with 26,527,622 to 51,801,664 reads per sample, and average quality >35.7 for all samples. After removing low-abundance transcripts, 357,203,972 reads (76.1%) mapped to 31,869 unique transcripts in the reference transcriptome (see Supplementary Table S2 for sample-specific read data). We assumed a negative binomial distribution for the count data and the log link function [77]. We used likelihood ratios based on Type III estimable functions to evaluate the significance of fixed effects. To adjust for multiple tests, we used the adaptive falsediscovery rates of Benjamini & Hochberg [78], as implemented in SAS Proc Multtest, and FDR<0.05 as the criterion for significance, except as noted. General and generalized linear model were implemented in SAS v. 9.3 [60] running under Linux 2.6.32._
