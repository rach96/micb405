## Methods

750 words

Reference your bash script(s) containing all commands. These should be available on a group member's 
GitHub page and/or emailed to Kim (kadm@mail.ubc.ca).
Explain why each step was used, what the outputs are, etc. 
For example, we began with paired-end reads and therefore, these 300 bp reads were combined into 
~300-350 bp contigs using the overlapping regions to yield a single contig for each sequence.
You need only explain the general flow of your pipeline since it is your results and discussion 
that will get into the meat of your specific parameter choices.


## Scripts for something...
script /dev/null
bash /home/micb405/Group12/Project3/scripts/len_400_t5_k8_diff2.sh

script /dev/null
bash /home/micb405/Group12/Project3/scripts/len_298_t5_k8_diff1.sh

script /dev/null
bash /home/micb405/Group12/Project3/scripts/distance_matrix_len_400_diff_1.sh


## FASTQC on raw
Change absolute path*

fastqc SI072_S3_150_1.fastq
fastqc SI072_S3_150_2.fastq

## Cutadapt: Remove 50 last bases
cutadapt -u -50 -o SI072_S3_150_1_trimmed.fastq SI072_S3_150_1.fastq
cutadapt -u -50 -o SI072_S3_150_2_trimmed.fastq SI072_S3_150_2.fastq

## FASTQC on trimmed
Change absolute path*

fastqc SI072_S3_150_1_trimmed.fastq
fastqc SI072_S3_150_2_trimmed.fastq

# Sequence Cleanup Steps:

## Create .files
make.file(inputdir=/home/micb405/Group12/Project3/fastq_files, prefix=Saanich150m)

## Make contigs
make.contigs(file=Saanich150m.files)

## Summarize output
summary.seqs(fasta=/home/micb405/Group12/Project3/contigs/Saanich150m.trim.contigs.fasta)

## Quality control
first:
screen.seqs(fasta=/home/micb405/Group12/Project3/contigs_2/Saanich150m.trim.contigs.fasta, group=Saanich150m.contigs.groups, maxambig=0, maxhomop=8, minlength=282, maxlength=298)

## Summarizing Resulting Sequences:
summary.seqs(fasta=Saanich150m.trim.contigs.good.fasta)

## Remove Duplicate Sequences:
set.dir(output=/home/micb405/Group12/Project3/screen_seqs/unique_seqs)
unique.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/Saanich150m.trim.contigs.good.fasta)

## Combine .names and.groups files into a count_table to shorten later steps that require both files
set.dir(output=/home/micb405/Group12/Project3/screen_seqs/unique_seqs)
count.seqs(name=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.names,group=/home/micb405/Group12/Project3/screen_seqs/Saanich150m.contigs.good.groups)

##Summarize again:
This summary will inform how well your quality control in screen.seqs went. Were you too strict? Do you have enough data left for coverage of your sample? Could you be more strict and get better quality data?

set.dir(output=/home/micb405/Group12/Project3/screen_seqs/unique_seqs)
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.fasta,count=Saanich150m.trim.contigs.good.count_table)


## Align sequences to a database. 
Default threshold: 0.5, default kmer length: 8

Miguel
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=5)

Rachel
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=8)

Amanda
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=12)

Kunye 
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.75, ksize=8)

## Summarize

Miguel (t0.5k5):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)

Rachel (t0.5k8):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)

Kunye (t0.75k8):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)

Amanda (t0.5k12)

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)


Some random commands: 

set.dir(output=/home/micb405/Group12/Project3/home/micb405/Group12/Project3/maxlength_400/align_seqs/de_replicated_post_cut)
summary.seqs(fasta=/home/micb405/Group12/Project3/maxlength_400/align_seqs/de_replicated_post_cut/Saanich150m.trim.contigs.good.unique.good.filter.unique.fasta, count=/home/micb405/Group12/Project3/maxlength_400/align_seqs/de_replicated_post_cut/Saanich150m.trim.contigs.good.unique.good.filter.count_table)


Cut the sequences to the same start and end. Use the previous summary to inform the cutoffs (stare and end) you want to use.

Miguel (t0.5k5):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)

Rachel (t0.5k8):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)

Amanda (t0.5k12):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)

Kunye (t0.75k8):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)

This summary will help you to see the result of the second screen.seqs step. If you've lost most of your sequences, then you need to alter the start/end parameters in the previous step.

Miguel (t0.5k5): 

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.good.count_table)

Rachel (t0.5k8):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/second_screen_output)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.good.count_table)

Amanda (t0.5k12):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.good.count_table)

Kunye (t0.75k8):

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.good.count_table)

## The alignment step adds . and - throughout the fasta file (hence why it is now a .align file instead of a .fasta). Columns of only . or - for all sequences provide no useful data and slow down later processes. So, we remove them here.
filter.seqs(fasta=, vertical=T, trump=.)

Miguel:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/editfasta)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

Rachel:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/edit_fasta_file)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

Amanda:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs/edit_fasta)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

Kunye:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/editfasta)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

## De-replicate again since alignment+filter may reveal additional identical sequences

Miguel:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/dereplicate)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/editfasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.good.count_table)

Rachel:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/edit_fasta_file/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.good.count_table)

Amanda:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12/dereplicate)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs/edit_fasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.good.count_table)

Kunye:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/dereplicate)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/editfasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.good.count_table)

## Pre-cluster very similar sequences. Consider the expected error rate of a 2x300 bp Illumina run. How many differences are actually expected sequencing error?

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/precluster)
pre.cluster(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.unique.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.count_table, diffs=1)

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2)
pre.cluster(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.unique.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.count_table, diffs=2)

Summarize Again

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

Pre-cluster Miguel t0.5k8 diff1:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

Pre-cluster Duncan t0.5k8 diff=2:

Pre-cluster: Rachel Diff 1

Pre-cluster: Rachel Diff 2

Diff1 chimera:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera)

chimera.vsearch(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

remove.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera)

summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

Diff 2 chimera:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera)

chimera.vsearch(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

remove.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

Rachel Diff1 chimera:

Rachel Diff2 chimera:

diff1:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton)
split.abund(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=1)

Diff2:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton)
split.abund(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=1)

## How many sequences were singletons?

Miguel diff1:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.abund.count_table)

Miguel Diff2:

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.abund.count_table)

Rachel Diff 1:

Rachel Diff 2:


Copy files over
system(cp Saanich.long.file.name.fasta Saanich.10m.final.fasta)

---------------------------------------------------------------------
# Clustering OTU:

## Calculate a distance matrix. Please use the lower triange (lt) format to save space.

400Diff1:
set.dir(output=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1)
dist.seqs(fasta=Saanich150m_t0.5k8_len400_diff1.fasta, output=lt)

298Diff1:

set.dir(output=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_1)
dist.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.fasta, output=lt)

298Diff2:

set.dir(output=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_2)
dist.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, output=lt)

## Cluster sequences. Carefully consider the method you want to use.  

set.dir(output=/home/micb405/Group12/Project3/clustering)

**400 Diff 1: 0.03**
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1/Saanich150m_t0.5k8_len400_diff1.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=opti)

**298 Diff 1: 0.03**
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_1/Saanich150m_t0.5k8diff1.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, method=opti)

**298 Diff 2: 0.03**
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_2/Saanich150m_t0.5k8diff2.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, method=opti)

Output File Names: 
/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list
/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.steps
/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.sensspec

**400 Diff 1: 0.01**
set.dir(output=/home/micb405/Group12/Project3/clustering/0.01)

cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1/Saanich150m_t0.5k8_len400_diff1.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=opti, cutoff=0.01)

**400 Diff 1: DGC**
set.dir(output=/home/micb405/Group12/Project3/clustering/dgc)

cluster(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=dgc, cutoff=0.03)


## Tutorial Command: Make Shared

set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files)

make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, label=0.03)

For Len298 Diff 1:

set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files)

make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, label=0.03)

For Len298 Diff 2:

set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files)

make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, label=0.03)

set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files/dgc)
make.shared(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, label=0.03)

Classify the sequences based on a database. Silva or GreenGenes? What bootstrap cutoff for confidence in your taxonomic assignment?

298diff2:

set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)

Doing: 400diff1

set.dir(output=/home/micb405/Group12/Project3/classifyseqs)

298diff1

set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)

classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)

400 diff1:

set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)

classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)

400diff2:

set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)

## Use the sequence classification output to apply the resulting taxonomy to OTUs

298diff2 (80-80):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, cutoff=80, threshold=80, basis=otu) 

298diff2 (60-60):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, cutoff=60, threshold=60, basis=otu) 

298diff1 (80-80):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, cutoff=80, threshold=80, basis=otu) 

298diff1 (60-60):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, cutoff=60, threshold=60, basis=otu) 

400 diff 1 (80-80):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 1 (60-60):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 

400 diff 2 (80-80):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 2 (60-60):

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)

classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, cutoff=60, threshold=60, basis=otu) 

400 diff 1 (80-80)0.01: 

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80/0.01)

classify.otu(list=/home/micb405/Group12/Project3/clustering/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 1 (60-60)0.01:

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60/0.01)

classify.otu(list=/home/micb405/Group12/Project3/clustering/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 

400 diff 1 (80-80)dgc:

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80/dgc)

classify.otu(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 1 (60-60)dgc:

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60/dgc)

classify.otu(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 


Summarize OTU table: 

298diff1:

set.dir(output=/home/micb405/Group12/Project3/summaryotu)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

298diff2:

set.dir(output=/home/micb405/Group12/Project3/summaryotu)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff2.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

400 diff 1:

set.dir(output=/home/micb405/Group12/Project3/summaryotu)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

400diff2:

set.dir(output=/home/micb405/Group12/Project3/summaryotu)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

400 diff 1 0.01:

set.dir(output=/home/micb405/Group12/Project3/summaryotu/0.01)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

set.dir(output=/home/micb405/Group12/Project3/summaryotu/0.01)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

400 diff 1 dgc:

set.dir(output=/home/micb405/Group12/Project3/summaryotu/dgc)

summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.shared, calc=nseqs-sobs-coverage-shannon-chao)

Python Script Commands:








