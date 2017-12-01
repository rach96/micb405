## Methods

750 words

Reference your bash script(s) containing all commands. These should be available on a group member's 
GitHub page and/or emailed to Kim (kadm@mail.ubc.ca).
Explain why each step was used, what the outputs are, etc. 
For example, we began with paired-end reads and therefore, these 300 bp reads were combined into 
~300-350 bp contigs using the overlapping regions to yield a single contig for each sequence.
You need only explain the general flow of your pipeline since it is your results and discussion 
that will get into the meat of your specific parameter choices.

## FASTQC on raw

```
fastqc /home/micb405/Group12/Project3/fastq_files/SI072_S3_150_1.fastq /home/micb405/Group12/Project3/fastqc
fastqc /home/micb405/Group12/Project3/fastq_files/SI072_S3_150_2.fastq /home/micb405/Group12/Project3/fastqc
```

## Cutadapt: Remove 50 last bases

```
cutadapt -u -50 -o SI072_S3_150_1_trimmed.fastq SI072_S3_150_1.fastq
cutadapt -u -50 -o SI072_S3_150_2_trimmed.fastq SI072_S3_150_2.fastq
```

## FASTQC on trimmed
Change absolute path*

```
fastqc /home/micb405/Group12/Project3/fastq_files/SI072_S3_150_1_trimmed.fastq /home/micb405/Group12/Project3/fastqc
fastqc /home/micb405/Group12/Project3/fastq_files/SI072_S3_150_2_trimmed.fastq /home/micb405/Group12/Project3/fastqc
```

# Sequence Cleanup Steps:

## Create .files

```
make.file(inputdir=/home/micb405/Group12/Project3/fastq_files, prefix=Saanich150m)
```

## Make contigs

```
make.contigs(file=Saanich150m.files)
```

## Summarize output

```
summary.seqs(fasta=/home/micb405/Group12/Project3/contigs/Saanich150m.trim.contigs.fasta)
```

## Quality control

```
screen.seqs(fasta=/home/micb405/Group12/Project3/contigs_2/Saanich150m.trim.contigs.fasta, group=Saanich150m.contigs.groups, maxambig=0, 
maxhomop=8, minlength=282, maxlength=298)
```

## Summarizing Resulting Sequences:

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.fasta)
```

## Remove Duplicate Sequences:

```
unique.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/Saanich150m.trim.contigs.good.fasta)
```

## Combine .names and.groups files into a count_table to shorten later steps that require both files

```
count.seqs(name=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.names,group=/home/micb405/Group12/Pro
ject3/screen_seqs/Saanich150m.contigs.good.groups)
```

## Summarize again:
This summary will inform how well your quality control in screen.seqs went. Were you too strict? Do you have enough data left for coverage 
of your sample? Could you be more strict and get better quality data?

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.fasta,count=Saanich150m.trim.contigs.good.count_table)
```
Note: The previous commands were executed for maxlengths of both 298 and 400

## Align sequences to a database. 
Default threshold: 0.5, default kmer length: 8
T0.5 represents a threshold of 0.5
k5 represents a kmer length of 5

**Miguel (t0.5k5):**

```
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, 
reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=5)
```

**Rachel (t0.5k8):**

```

align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, 
reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=8)
```

**Amanda (t0.5k12):**

```
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, 
reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=12)
```

**Kunye (t0.75k8):**

```
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, 
reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.75, ksize=8)
```

## Summarize

**Miguel (t0.5k5):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)
```

**Rachel (t0.5k8):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)
```

**Amanda (t0.5k12):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)
```

**Kunye (t0.75k8):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)
```


## Keep Sequences With Same Start and End Alignment Positions

**Miguel (t0.5k5):**

```
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, 
summary=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
```

**Rachel (t0.5k8):**

```
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, 
summary=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
```

**Amanda (t0.5k12):**

```
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, 
summary=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
```

**Kunye (t0.75k8):**

```
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.align, 
count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, 
summary=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
```

## Summarize Previous Step

**Miguel (t0.5k5):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.good.count_table)
```

**Rachel (t0.5k8):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.good.align, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.good.count_table)
```

**Amanda (t0.5k12):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.good.align, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.good.count_table)
```

**Kunye (t0.75k8):**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, 
count=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.good.count_table)
```

## Change .align File to .fasta File

**Miguel (t0.5k5):**

```
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, 
trump=.)
```

**Rachel (t0.5k8):**

```
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```

**Amanda (t0.5k12):**

```
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```

**Kunye (t0.75k8):**

```
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, 
trump=.)
```

## Remove Duplicate Sequences

**Miguel (t0.5k5):**

```
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/editfasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta
,count=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.good.count_table)
```

**Rachel (t0.5k8):**

```
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/edit_fasta_file/Saanich150m.trim.contigs.good.unique.good.filter.fasta,c
ount=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.good.count_table)
```

**Amanda (t0.5k12):**

```
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs/edit_fasta/Saanich150m.trim.contigs.good.unique.good.filter.fas
ta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.good.count_table)
```

**Kunye (t0.75k8):**

```
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/editfasta/Saanich150m.trim.contigs.good.unique.good.filter.fasa
,count=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.good.count_table)
```

## Pre-cluster Sequences With 1 and 2 Mismatches

**298 Diff 1**

```
pre.cluster(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.uni
que.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.count_table, 
diffs=1)
```

```
pre.cluster(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.uni
que.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again/Saanich150m.trim.contigs.good.unique.good.filter.count_table, 
diffs=2)
```

## Summarize The Previous Step

## Pre-Cluster Summary:

**Pre-cluster Miguel t0.5k8 diff1:**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique
.precluster.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.c
ount_table)
```

**Pre-cluster Duncan t0.5k8 diff=2:**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique
.precluster.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.c
ount_table)
```


**Pre-cluster: Rachel Diff 1**

```
pre.cluster(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.count_table, diffs=1)
```

**Pre-cluster: Rachel Diff 2**

```
pre.cluster(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.count_table, diffs=2)
```

## Pre-cluster Summary

**298 Diff1**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filte
r.unique.precluster.pick.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.prec
luster.denovo.vsearch.pick.count_table)
```

**298 Diff2**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filte
r.unique.precluster.pick.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.prec
luster.denovo.vsearch.pick.count_table)
```

**400 Diff1**

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
```

**400 Diff2**

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
```


## Removing Chimeras

**Diff1 chimera:**

```
chimera.vsearch(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.uni
que.precluster.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.c
ount_table, dereplicate=t)
```

```
remove.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/Saanich150m.trim.contigs.good.unique.good.filter.unique.
precluster.fasta, 
accnos=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.pre
cluster.denovo.vsearch.accnos)
```


**Diff 2 chimera:**

```
chimera.vsearch(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.uni
que.precluster.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.c
ount_table, dereplicate=t)
```

```
remove.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/Saanich150m.trim.contigs.good.unique.good.filter.unique.
precluster.fasta, 
accnos=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.pre
cluster.denovo.vsearch.accnos)
```

**400 Diff1 chimera:**

```
chimera.vsearch(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
```

```
remove.seqs(fasta=/home/micb405/Group12/Project3/maxlength_400/align_seqs/pre_clustering/diff_1/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```


**400 Diff2 chimera:**

```
chimera.vsearch(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
```

```
remove.seqs(fasta=/home/micb405/Group12/Project3/maxlength_400/align_seqs/pre_clustering/diff_2/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```

## Chimera Summary:

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
```
* Note: This command was run for all 4 conditions (298 diff1/2 and 400 diff1/2)


## Removing Singletons

**298 diff1:**

```
split.abund(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter
.unique.precluster.pick.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.prec
luster.denovo.vsearch.pick.count_table, cutoff=1)
```

**298 Diff2:**

```
split.abund(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter
.unique.precluster.pick.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.prec
luster.denovo.vsearch.pick.count_table, cutoff=1)
```

**400 Diff1:**

```
split.abund(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=1)
```


**400 Diff2:**

```
split.abund(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=1)
```


## Singleton Summary?

**298 diff1:**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton/Saanich150m.trim.contigs.good.unique.good.fil
ter.unique.precluster.pick.abund.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton/Saanich150m.trim.contigs.good.unique.good.filter.unique.pr
ecluster.denovo.vsearch.pick.abund.count_table)
```

**298 Diff2:**

```
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton/Saanich150m.trim.contigs.good.unique.good.fil
ter.unique.precluster.pick.abund.fasta, 
count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton/Saanich150m.trim.contigs.good.unique.good.filter.unique.pr
ecluster.denovo.vsearch.pick.abund.count_table)
```

**400 Diff 1:**

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.abund.count_table)
```

**400 Diff 2:**

```
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta, count=Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.abund.count_table)
```


---------------------------------------------------------------------
# Clustering OTU:

## Calculate a distance matrix. Please use the lower triange (lt) format to save space.

**400Diff1:**

```
dist.seqs(fasta=Saanich150m_t0.5k8_len400_diff1.fasta, output=lt)
```

**298Diff1:**

```
dist.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.fasta, output=lt)
```

**298Diff2:**

```
dist.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, output=lt)
```

## Cluster sequences. Carefully consider the method you want to use.  

**298 Diff 1: 0.03**
```
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_1/Saanich150m_t0.5k8diff1.phylip.dist, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, method=opti)
```

**298 Diff 2: 0.03**
```
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_2/Saanich150m_t0.5k8diff2.phylip.dist, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, method=opti)
```

**400 Diff 1: 0.03**
```
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1/Saanich150m_t0.5k8_len400_diff1.phylip.dist, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=opti)
```

**400 Diff 2: 0.03**
```
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_2/Saanich150m_t0.5k8_len400_diff2.phylip.dist, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, method=opti)
```

**400 Diff 2: 0.03**

```
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_2/Saanich150m_t0.5k8_len400_diff2.phylip.dist, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, method=opti, cutoff=0.01)
```

**400 Diff 2: DGC**

```
cluster(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.fasta, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, method=dgc, cutoff=0.03)
```

## Tutorial Command: Make Shared

```
make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, label=0.03)
```

**For Len 298 Diff 1:**

```
make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, label=0.03)
```

**For Len 298 Diff 2:**

```
make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, label=0.03)
```

**For Len 400 Diff 2:**

```
make.shared(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, label=0.03)
```

---------------------------------------------------------------------

## Classify sequences based on a Silva

**298 diff 1:**

```
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.fasta, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, 
template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, 
cutoff=60)
```

**298 diff 2:**

```
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, 
template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, 
cutoff=60)
```

**400 diff 1:**

```
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, 
template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, 
cutoff=60)
```

**400 diff 2:**

```
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, 
template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, 
cutoff=60)
```

## Use the sequence classification output to apply the resulting taxonomy to OTUs: Label 0.03

**298 diff 1 (60-60):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, cutoff=60, threshold=60, basis=otu) 
```

**298 diff 1 (80-80):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, cutoff=80, threshold=80, basis=otu) 
```

**298 diff 2 (80-80):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff2.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, cutoff=80, threshold=80, basis=otu) 
```

**298 diff 2 (60-60):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff2.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, cutoff=60, threshold=60, basis=otu) 
```

**400 diff 1 (60-60):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 
```

**400 diff 1 (80-80):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 
```

**400 diff 2 (60-60):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, cutoff=60, threshold=60, basis=otu) 
```

**400 diff 2 (80-80):**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, cutoff=80, threshold=80, basis=otu)
```

## Use the sequence classification output to apply the resulting taxonomy to OTUs: Label 0.01

**400 diff 2 (60-60)0.01:**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu)
```

**400 diff 2 (80-80)0.01:**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 
```

## Use the sequence classification output to apply the resulting taxonomy to OTUs: Clustering DGC

**400 diff 2 (60-60)dgc:**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 
```

**400 diff 2 (80-80)dgc:**

```
classify.otu(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, 
taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, 
count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 
```

## Summarize OTU table: 

**298 diff 1:**

```
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff1.phylip.opti_mcc.shared, calc=nseqs-
sobs-coverage-shannon-chao)
```

**298 diff 2:**

```
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff2.phylip.opti_mcc.shared, calc=nseqs-
sobs-coverage-shannon-chao)
```

**400 diff 1:**

```
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, 
calc=nseqs-sobs-coverage-shannon-chao)
```

**400 diff 2:**

```
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.shared, 
calc=nseqs-sobs-coverage-shannon-chao)
```

**400 diff 2 0.01:**

```
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/0.01/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.shared, 
calc=nseqs-sobs-coverage-shannon-chao)
```

**400 diff 2 dgc:**

```
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/dgc/Saanich150m_t0.5k8_len400_diff2.dgc.unique_list.shared, 
calc=nseqs-sobs-coverage-shannon-chao)
```

## Python Script Command:

```
python taxonomy1b.py
```

* Note: these commands were executed on all of the 400 diff 2 files








