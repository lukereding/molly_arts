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

### create a job file
`touch filter_job`<br>`for i in *.fastq; do echo "fastq_quality_filter -q 20 -p 80 -i $i -Q 33 -o $i.filtered" >> filter_job; done`<br>`cat filter_job`       

>fastq_quality_filter -q 20 -p 80 -i SRR1161450_1.fastq -Q 33 -o SRR1161450_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1161451_1.fastq -Q 33 -o SRR1161451_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1165201_1.fastq -Q 33 -o SRR1165201_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1165203_1.fastq -Q 33 -o SRR1165203_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166366_1.fastq -Q 33 -o SRR1166366_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166367_1.fastq -Q 33 -o SRR1166367_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166368_1.fastq -Q 33 -o SRR1166368_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166369_1.fastq -Q 33 -o SRR1166369_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166370_1.fastq -Q 33 -o SRR1166370_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166371_1.fastq -Q 33 -o SRR1166371_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166372_1.fastq -Q 33 -o SRR1166372_1.fastq.filtered

That should work fine as a job file.

`module load fastx_toolkit` `module load launcher`       

### submit the job
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

### getting a transcriptome

I first decided to use the guppy transcriptome:       
`wget ftp://ftp.tuebingen.mpg.de/ebio/publication_data/esharma/guppy_trans/trin_cuff_v14_cdhit90.fa.gz`   # download the transcriptome     
`gunzip trin_cuff_v14_cdhit90.fa.gz` #unzip it     
`mv trin_cuff_v14_cdhit90.fa guppy_transcriptome` # rename it    


### making transcriptome annotations
#### doing the blasts
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


--------------------
#### processing blast results
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
## mapping

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
## SAM conversion, sorting, and indexing

We now convert the SAM files that resulted from the previous step to BAM files, which are less unwieldy. 

`for file in SRR*.sam; do echo "samtools view -b -S $file > ${file%.sam}.bam && samtools sort ${file%.sam}.bam ${file%.sam}.bam.sort && samtools index ${file%.sam}.bam" >> bam_commands; done`

Launch it:     
`launcher_creator.py -j bam_commands -n bam_commands_job -l bam_commands_job -a Sailfin_RNASeq -e lukereding@utexas.edu -q normal -t 2:00:00` 
`qsub bam_commands_job`


Now we can assess mapping quality for the reads from each individual like:

`samtools flagstat SRR1165201_1.bam`

-----------------

## gene counting


`module load samtools`        

For each individual (i.e. each BAM file), we can use idxstats in `samtools` to extract a list of genes and their counts:
`samtools idxstats SRR1166372_1.bam.sort.bam | cut -f 1,3`

> **note: check this.** initial tests showed that this shows some duplicate genes. see the bsa hscript [here](https://www.biostars.org/p/14531/) for a possible, better alternative.


---------------

### from the supplemental materials:
_*for reference:*_         
_(a) RNA-seq data collection and analysis RNA was extracted from whole brain tissue using RNeasy Lipid Tissue Mini Kit (Qiagen Hilden, Germany). We prepared separate RNA sequencing libraries from whole brains for each individual, using unique index sequences from the Illumina Tru-Seq RNA kit following manufacturers instructions. Sequencing libraries were constructed and sequenced on three lanes of an Illumina HiSeq 2000 at the HudsonAlpha Genomic Services Laboratory (Huntsville, Alabama) in April 2012. The reference P. reticulata assembly was constructed from a data set containing > 450 million 100-bp paired end reads, which were filtered for high quality sequence and normalized in-silico to compress the range in kmer abundance. We used SeqMan NGEN 4.1.2 (Madison, WI) [74,75] to perform the assembly. Contigs from the assembly were annotated by blastx queries against SwissProt (database downloaded Oct 6, 2012), UniProt/Trembl (Nov 28, 2012), and nr (Dec 11, 2012). Default parameters were used in the blastx queries, with e-value cutoff of 10-4 . The assembly and individual reads will be deposited in a publicly-available archive before publication. Reads were mapped to the reference assembly using Bowtie 2 v 2.0.0 on a server running Red Hat Enterprise Linux 6.1. We used a seed size of 20 bp, with no mismatches allowed in the seed (run options: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75). We retained mappings with quality scores > 30 (< 0.001 probability that the read maps elsewhere in the reference), and applied an abundance filter to retain only transcripts represented by more than 1 count per million reads in at least three samples. This filtering resulted in 31,869 transcripts remaining in the data set. We used the number of reads mapping to each of those transcripts, along with TMM-normalized library sizes [76] to analyse differential expression. We obtained 462,537,724 100-bp reads that passed the machine quality filter, with 26,527,622 to 51,801,664 reads per sample, and average quality >35.7 for all samples. After removing low-abundance transcripts, 357,203,972 reads (76.1%) mapped to 31,869 unique transcripts in the reference transcriptome (see Supplementary Table S2 for sample-specific read data). We assumed a negative binomial distribution for the count data and the log link function [77]. We used likelihood ratios based on Type III estimable functions to evaluate the significance of fixed effects. To adjust for multiple tests, we used the adaptive falsediscovery rates of Benjamini & Hochberg [78], as implemented in SAS Proc Multtest, and FDR<0.05 as the criterion for significance, except as noted. General and generalized linear model were implemented in SAS v. 9.3 [60] running under Linux 2.6.32._
