# playing around with sailfin RNA-seq data from Fraser et al. 2014
## getting the data
The data are stored on NCBI's short read archive, which is a terrible site. You cannot press a button or copy a link to download the raw data. I had to follow instructions [here](https://www.biostarhandbook.com/unit/sra/short-read-archive.html#sra). I first installed _sratoolkit_ on my Mac using Homebrew like `brew install sratoolkit`. I could then use the `prefetch` command to download the reads corresponding to each individual. I put the names of all the runs in a text document called runs, then I wrote a script, `download_short_read_in_par.sh`, to download the data in parrallel (this script and others are in this repository). Once the data are downloaded, they are in `.sra` format and not `.fastq` like you need them to be. So I wrote a second script, `sra_to_fastq.sh`, to convert each `.sra` file to `.fastq`, which takes a bit of time. Finally, I transferred all 11 fastq files to lonestar using `scp /path/to/fastq/files/*.fastq reding@lonestar.tacc.utexas.edu`.         

The reads are 100 bp single end reads. `head -4 SRR1161450_1.fastq` gives:

@SRR1161450.1 HWI-ST619:197:D0VKFACXX:6:1101:1073:1995 length=100 CAGNTTTAGTCCAAAGTTTCTATATACAGTCAGAGATGAAACAGTTCTGGGCTTGGCCAAGCTGAAAAGAGGCTTCAGCTCCAGCTGAGTTCATCATTTN +SRR1161450.1 HWI-ST619:197:D0VKFACXX:6:1101:1073:1995 length=100 B@C#4ADDHHHHHIJJHIJIJFIIJJJIJIJJJJJIJJHIJJJJGIJJJJJJJJJJJJJJJJJJHGIJJIJJHFHHFFFFFFEECEECB;CAACDDDEC#

--------------------------------------------------------------------------------

## quality filtering
In the paper, they don't say if or how they did quality filtering. Totally arbitrarily, I decided to filter the reads so that at least 80% of each read has a quality score of 30 or above.

`pwd`<br>/scratch/02535/reding/molly_arts     

### create a job file
`touch filter_job`<br>`for i in *.fastq; do echo "fastq_quality_filter -q 20 -p 80 -i $i -Q 33 -o $i.filtered" >> filter_job; done`<br>`cat filter_job`       

fastq_quality_filter -q 20 -p 80 -i SRR1161450_1.fastq -Q 33 -o SRR1161450_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1161451_1.fastq -Q 33 -o SRR1161451_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1165201_1.fastq -Q 33 -o SRR1165201_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1165203_1.fastq -Q 33 -o SRR1165203_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166366_1.fastq -Q 33 -o SRR1166366_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166367_1.fastq -Q 33 -o SRR1166367_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166368_1.fastq -Q 33 -o SRR1166368_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166369_1.fastq -Q 33 -o SRR1166369_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166370_1.fastq -Q 33 -o SRR1166370_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166371_1.fastq -Q 33 -o SRR1166371_1.fastq.filtered fastq_quality_filter -q 20 -p 80 -i SRR1166372_1.fastq -Q 33 -o SRR1166372_1.fastq.filtered

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
originally: 35832124 , filtered: 32671630 >> filtering_results         
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

The transcriptome is just a list of sequences that are likely genes; we don't know the identites of the genes. For this, we have to do a series of BLASTs, blasting each sequence against know gene sequences and names our genes based on the results. The Swiss-Prot database contains a bunch of sequences that we can blast against (I think). 


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


-------------------
#### (next steps--from Misha's pipeline)

unzipping      
gunzip uniprot_sprot.fasta.gz &         
 [1] 17275         
gunzip idmapping_selected.tab.gz &         
[2] 17278         
      
indexing the fasta database         
module load blast             
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb            
launcher_creator.py -j mdb -n mdb -l mmm            
qsub mmm           
 
splitting the transcriptome into 40 chunks     
splitFasta.pl monti_coral_iso.fasta 40        

blasting all 40 chunks to uniprot in parallel, 3 cores per chunk         
module load blast           
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta              -evalue 0\.0001 -num_threads 3 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl            
launcher_creator.py -j bl -n blast -l blj -q normal -t 24:00:00        
cat blj | perl -pe 's/12way .+$/4way 120/' >bljj           
qsub bljj           
 
watching progress:         
grep "Query= " subset*.br | wc -l               

if the blast did not finish in time, try splitting the transcriptome into 120 parts initially        

combining all blast results       
cat subset*br > myblast.br         
rm subset*br            

if you have no assembler-derived isogroups, use cd-hit-est to cluster contigs.        
to look for 99% or better matches between contigs taking 30% of length of either longer or shorter sequence:         
cd-hit-est -i monti_coral_iso.fasta -o transcriptome_clust.fasta -c 0.99 -G 0 -aL 0.3 -aS 0.3
-----------------

---------------

### from the supplemental materials:
_(a) RNA-seq data collection and analysis RNA was extracted from whole brain tissue using RNeasy Lipid Tissue Mini Kit (Qiagen Hilden, Germany). We prepared separate RNA sequencing libraries from whole brains for each individual, using unique index sequences from the Illumina Tru-Seq RNA kit following manufacturers instructions. Sequencing libraries were constructed and sequenced on three lanes of an Illumina HiSeq 2000 at the HudsonAlpha Genomic Services Laboratory (Huntsville, Alabama) in April 2012. The reference P. reticulata assembly was constructed from a data set containing > 450 million 100-bp paired end reads, which were filtered for high quality sequence and normalized in-silico to compress the range in kmer abundance. We used SeqMan NGEN 4.1.2 (Madison, WI) [74,75] to perform the assembly. Contigs from the assembly were annotated by blastx queries against SwissProt (database downloaded Oct 6, 2012), UniProt/Trembl (Nov 28, 2012), and nr (Dec 11, 2012). Default parameters were used in the blastx queries, with e-value cutoff of 10-4 . The assembly and individual reads will be deposited in a publicly-available archive before publication. Reads were mapped to the reference assembly using Bowtie 2 v 2.0.0 on a server running Red Hat Enterprise Linux 6.1. We used a seed size of 20 bp, with no mismatches allowed in the seed (run options: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75). We retained mappings with quality scores > 30 (< 0.001 probability that the read maps elsewhere in the reference), and applied an abundance filter to retain only transcripts represented by more than 1 count per million reads in at least three samples. This filtering resulted in 31,869 transcripts remaining in the data set. We used the number of reads mapping to each of those transcripts, along with TMM-normalized library sizes [76] to analyse differential expression. We obtained 462,537,724 100-bp reads that passed the machine quality filter, with 26,527,622 to 51,801,664 reads per sample, and average quality >35.7 for all samples. After removing low-abundance transcripts, 357,203,972 reads (76.1%) mapped to 31,869 unique transcripts in the reference transcriptome (see Supplementary Table S2 for sample-specific read data). We assumed a negative binomial distribution for the count data and the log link function [77]. We used likelihood ratios based on Type III estimable functions to evaluate the significance of fixed effects. To adjust for multiple tests, we used the adaptive falsediscovery rates of Benjamini & Hochberg [78], as implemented in SAS Proc Multtest, and FDR<0.05 as the criterion for significance, except as noted. General and generalized linear model were implemented in SAS v. 9.3 [60] running under Linux 2.6.32._
