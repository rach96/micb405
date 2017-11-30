Methods

750 words

Reference your bash script(s) containing all commands. These should be available on a group member's 
GitHub page and/or emailed to Kim (kadm@mail.ubc.ca).
Explain why each step was used, what the outputs are, etc. 
For example, we began with paired-end reads and therefore, these 300 bp reads were combined into 
~300-350 bp contigs using the overlapping regions to yield a single contig for each sequence.
You need only explain the general flow of your pipeline since it is your results and discussion 
that will get into the meat of your specific parameter choices.

Examples taken from: https://www.mothur.org/wiki/MiSeq_SOP

Steps To Examine: https://www.mothur.org/wiki/Category:Commands
Sequence Alignment (Step before also)
Chimera Identification: 
Cluster sequences: 

script /dev/null
bash /home/micb405/Group12/Project3/scripts/len_400_t5_k8_diff2.sh


script /dev/null
bash /home/micb405/Group12/Project3/scripts/len_298_t5_k8_diff1.sh

script /dev/null
bash /home/micb405/Group12/Project3/scripts/distance_matrix_len_400_diff_1.sh

FASTQC on raw
Change absolute path*

fastqc SI072_S3_150_1.fastq
fastqc SI072_S3_150_2.fastq

Cutadapt
http://cutadapt.readthedocs.io/en/stable/guide.html#modifying-reads

(From the folder with the fastq files or change absolute path)

Remove 50 last bases
cutadapt -u -50 -o SI072_S3_150_1_trimmed.fastq SI072_S3_150_1.fastq
cutadapt -u -50 -o SI072_S3_150_2_trimmed.fastq SI072_S3_150_2.fastq

FASTQC on trimmed
Change absolute path*

fastqc SI072_S3_150_1_trimmed.fastq
fastqc SI072_S3_150_2_trimmed.fastq

Sequence Cleanup Steps:

Create .files
make.file(inputdir=/home/micb405/Group12/Project3/fastq_files, prefix=Saanich150m)

Make contigs
make.contigs(file=Saanich150m.files)

Summarize output
summary.seqs(fasta=/home/micb405/Group12/Project3/contigs/Saanich150m.trim.contigs.fasta)


Quality control
What cutoffs do you use given the 2x300 bp quality? How long do we expect the 16S V4-V5 region to be?



first:
screen.seqs(fasta=/home/micb405/Group12/Project3/contigs_2/Saanich150m.trim.contigs.fasta, group=Saanich150m.contigs.groups, maxambig=0, maxhomop=8, minlength=282, maxlength=298)

NOTE: Output is in screen_seqs folder
Maxambig = 0?
While we don't necessarily know the longest acceptable homopolymer for a 16S rRNA gene, the max length of 31 is clearly a sequencing artifact. If you are interested in removing sequences with excessively long homopolymers, then you should use the maxhomop option:
Maxhomop = 8 ?
Minlength = 295 ?
Maxlength = 300 ?
Group = Saanich150m.groups ?

Summarizing Resulting Sequences:
summary.seqs(fasta=Saanich150m.trim.contigs.good.fasta)



400_len_diff_1




Remove Duplicate Sequences:

My code:
set.dir(output=/home/micb405/Group12/Project3/screen_seqs/unique_seqs)
unique.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/Saanich150m.trim.contigs.good.fasta)

Combine .names and.groups files into a count_table to shorten later steps that require both files

My Code:
set.dir(output=/home/micb405/Group12/Project3/screen_seqs/unique_seqs)
count.seqs(name=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.names,group=/home/micb405/Group12/Project3/screen_seqs/Saanich150m.contigs.good.groups)


Summarize again. This summary will inform how well your quality control in screen.seqs went. Were you too strict? Do you have enough data left for coverage of your sample? Could you be more strict and get better quality data?


400_len_1_diff_1







My Code:
set.dir(output=/home/micb405/Group12/Project3/screen_seqs/unique_seqs)
summary.seqs(fasta=Saanich150m.trim.contigs.good.unique.fasta,count=Saanich150m.trim.contigs.good.count_table)



TODO:

Align sequences to a database. Please only use Silva in this step as this is a LONG process and alignments to GreenGenes are of poorer quality. We recommend using flip=T. What does this parameter mean and why might we want to include it?
Nohup equivalent: input script /dev/null and then run your shell script bash my_script.sh
Default threshold: 0.5, default kmer length: 8

Fasta:/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta


Miguel
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=5)

Rachel
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=8)

mothur "#set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8); ./t5_k8.sh"

Amanda
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.50, ksize=12)

Kunye 
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8)
align.seqs(fasta=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T, threshold=0.75, ksize=8)

Change kmer size?? 5, 8 & 12 -> fastest alignment is usually the best (CHECK LOGS!)
Threshold: do we want to do 0.5 and 0.75? (let’s start with 0.5 since it’s default)




Options:

flip and threshold
The threshold parameter is used to specify a cutoff at which an alignment is deemed 'bad' and the reverse complement may be tried. The default threshold is 0.50, meaning if 50% of the bases are removed in the alignment process. 
The flip parameter is used to specify whether or not you want mothur to try the reverse complement of a sequence if the sequence falls below the threshold. The default is false. If the flip parameter is set to true the reverse complement of the sequence is aligned and the better alignment is reported.
In the pairwise alignment portion of the aligning procedure, the default reward for a match is +1 and the penalties for a mismatch, opening and extending a gap are -1, -2, and -1. Our experience has shown that these produce the best alignments for 16S rRNA gene sequences. 
More parameters at: https://www.mothur.org/wiki/Align.seqs





Summarize


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

overview
set.dir(output=/home/micb405/Group12/Project3/alignseqs/$PATH)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/$PATH, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table)
Count path:

set.dir(output=/home/micb405/Group12/Project3/home/micb405/Group12/Project3/maxlength_400/align_seqs/de_replicated_post_cut)
summary.seqs(fasta=/home/micb405/Group12/Project3/maxlength_400/align_seqs/de_replicated_post_cut/Saanich150m.trim.contigs.good.unique.good.filter.unique.fasta, count=/home/micb405/Group12/Project3/maxlength_400/align_seqs/de_replicated_post_cut/Saanich150m.trim.contigs.good.unique.good.filter.count_table)

Count path:


Cut the sequences to the same start and end. Use the previous summary to inform the cutoffs (stare and end) you want to use.


screen.seqs(fasta=, count=, summary=, start=, end=)
Miguel (t0.5k5)
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.5k5/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
Rachel (t0.5k8):
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
Amanda (t0.5k12)
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
Kunye (t0.75k8)
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs)
screen.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.align, count=/home/micb405/Group12/Project3/screen_seqs/unique_seqs/Saanich150m.trim.contigs.good.count_table, summary=/home/micb405/Group12/Project3/alignseqs/t0.75k8/Saanich150m.trim.contigs.good.unique.summary, start=10370, end=25318)
This summary will help you to see the result of the second screen.seqs step. If you've lost most of your sequences, then you need to alter the start/end parameters in the previous step.
Miguel (t0.5k5): 

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.good.count_table)
Rachel (t0.5k8) :

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/second_screen_output)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.good.count_table)

Amanda

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.good.count_table)


Kunye (t0.75k8)

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs)
summary.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, count=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.good.count_table)

The alignment step adds . and - throughout the fasta file (hence why it is now a .align file instead of a .fasta). Columns of only . or - for all sequences provide no useful data and slow down later processes. So, we remove them here.
filter.seqs(fasta=, vertical=T, trump=.)


Additional info: Now we know our sequences overlap the same alignment coordinates, we want to make sure they only overlap that region. So we'll filter the sequences to remove the overhangs at both ends. Since we've done paired-end sequencing, this shouldn't be much of an issue, but whatever. 
This means that our initial alignment was 13425 columns wide and that we were able to remove 13049 terminal gap characters using trump=. and vertical gap characters using vertical=T. The final alignment length is 376 columns. Because we've perhaps created some redundancy across our sequences by trimming the ends, we can re-run unique.seqs:

Miguel:
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/editfasta)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

Rachel:
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/edit_fasta_file)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

Amanda
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs/edit_fasta)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

Kunye
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/editfasta)
filter.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.unique.good.align, vertical=T, trump=.)

De-replicate again since alignment+filter may reveal additional identical sequences



unique.seqs(fasta=, count=)
Miguel:
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k5/dereplicate)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/editfasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k5/cutseqs/Saanich150m.trim.contigs.good.good.count_table)


Rachel:
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/de_replicated_again)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/edit_fasta_file/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/Saanich150m.trim.contigs.good.good.count_table)

Amanda
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k12/dereplicate)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k12/cutseqs/edit_fasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.5k12/Saanich150m.trim.contigs.good.good.count_table)

Kunye

set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.75k8/dereplicate)
unique.seqs(fasta=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/editfasta/Saanich150m.trim.contigs.good.unique.good.filter.fasta,count=/home/micb405/Group12/Project3/alignseqs/t0.75k8/cutseqs/Saanich150m.trim.contigs.good.good.count_table)

Pre-cluster very similar sequences. Consider the expected error rate of a 2x300 bp Illumina run. How many differences are actually expected sequencing error?
Pitfalls of clustering: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909393/

Error rate for Illumina Miseq = ~0.1%
Primary errors: substitution
“For Illumina, errors at ≤0.1% is achieved for ≥ 75-85% of bases. In general, Illumina keeps similar error rate criteria, but extends the maximum read length with new versions of their chemistry.”
Source: http://onlinelibrary.wiley.com.ezproxy.library.ubc.ca/doi/10.1111/j.1755-0998.2011.03024.x/full 
http://www.molecularecologist.com/next-gen-table-3c-2014/ 

Two examples of diffs = 2
http://www.microbiota.org/cgi-bin/tmp/mothur.cgi 
Kim’s work: https://rpubs.com/dillmcfarlan/mothurSOP 

From the readings: https://github.com/EDUCE-UBC/MICB405_project3/blob/master/Readings/Loman%20NJ%202012%20Nature%20Biotech_seq%20platforms.pdf 
“MiSeq produced the highest quality reads, owing to a low substitution error rate (0.1 substitutions per 100 bases)
Indel = insertion or deletion of bases in the genome
“Indels were detected very infrequently in MiSeq data with <0.001 indels per 100 bases.”

Notes: The next thing we want to do to further de-noise our sequences is to pre-cluster the sequences using the pre.cluster command allowing for up to 2 differences between sequences. 
This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged. We generally favor allowing 1 difference for every 100 bp of sequence:
Diffs (Parameter)
By default the pre.cluster command will look for sequences that are within 1 mismatch of the sequence being considered. With the diffs option you can change this threshold. For example:


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



Vsearch Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5075697/
Chimera Identification and Removal:
Tutorial: chimera.vsearch
We now have 5721 unique sequences. At this point we have removed as much sequencing error as we can and it is time to turn our attention to removing chimeras. We'll do this using the VSEARCH algorithm that is called within mothur using the chimera.vsearch command. Again, this command will split the data by sample and check for chimeras. Our preferred way of doing this is to use the abundant sequences as our reference. In addition, if a sequence is flagged as chimeric in one sample, the default (dereplicate=F) is to remove it from all samples. Our experience suggests that this is a bit aggressive since we've seen rare sequences get flagged as chimeric when they're the most abundant sequence in another sample. This is how we do it:




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

Using stability.trim.contigs.good.good.count_table as input file for the count parameter.
Using stability.trim.contigs.good.unique.good.align as input file for the fasta parameter.
Notes: Note that we went from 128,655 to 118,150 sequences for a reduction of 8.1%; this is a reasonable number of sequences to be flagged as chimeric. 
To facilitate downstream analyses, we will remove sequences that occur only once in the whole dataset (e.g.singletons). We need to share the server and decrease each group's usage.
diff1
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/singleton)
split.abund(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff1/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=1)
Diff2
set.dir(output=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/singleton)
split.abund(fasta=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=/home/micb405/Group12/Project3/alignseqs/t0.5k8/preclusterdiff2/chimera/Saanich150m.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=1)


How many sequences were singletons?
summary.seqs(fasta=, count=)
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
At this point the tutorial diverges from the steps in the workflow so I’m not sure if any of them are correct…. 
Clustering OTU’s Steps: Moving forward, you should use your final files.
Tutorial:


---------------------------------------------------------------------
Calculate a distance matrix. Please use the lower triange (lt) format to save space.
set.dir(output=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1)
dist.seqs(fasta=Saanich150m_t0.5k8_len400_diff1.fasta, output=lt)
298Diff1:
set.dir(output=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_1)
dist.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.fasta, output=lt)
298Diff2:
set.dir(output=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_2)
dist.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, output=lt)

Cluster sequences. Carefully consider the method you want to use.  
Cutoff? 
set.dir(output=/home/micb405/Group12/Project3/clustering)
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1/Saanich150m_t0.5k8_len400_diff1.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=opti)
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_1/Saanich150m_t0.5k8diff1.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, method=opti)
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_298_diff_2/Saanich150m_t0.5k8diff2.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, method=opti)

Output File Names: 
/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list
/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.steps
/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.sensspec

set.dir(output=/home/micb405/Group12/Project3/clustering/0.01)
cluster(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1/Saanich150m_t0.5k8_len400_diff1.phylip.dist, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=opti, cutoff=0.01)
set.dir(output=/home/micb405/Group12/Project3/clustering/dgc)
cluster(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, method=dgc, cutoff=0.03)
Cutoff (https://mothur.org/wiki/Cluster#cutoff)
“With the opticlust method the list file is created for the cutoff you set. The default cutoff is 0.03. With the average neighbor, furthest neighbor and nearest neighbor methods the cutoff should be significantly higher than the desired distance in the list file. We suggest cutoff=0.20. This will provide a boost in speed and less RAM will be required than if you didn't set the cutoff for reading in the matrix. The cutoff can be set for the cluster command as follows:”
cluster(phylip=, count=, method=dgc)  (if we have time)
http://msphere.asm.org/content/2/2/e00073-17.figures-only
Use OptiClust and DGC VSEARCH (if time)


Options: The methods available in mothur include opticlust (opti), average neighbor (average), furthest neighbor (furthest), nearest neighbor (nearest), Vsearch agc (agc), Vsearch dgc (dgc). By default cluster() uses the opticlust algorithm; this can be changed with the method option.
Default method: opticlust: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5343174/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5343174/figure/fig3/
Consistently, the average neighbor algorithm was identified as among the best or as the best algorithm. Other hierarchical algorithms, such as furthest and nearest neighbor, which do not permit the formation of FPs or FNs, respectively, fared significantly worse. The distance-based greedy clustering as implemented in VSEARCH has also performed well. The computational resources required to complete the average neighbor algorithm can be significant for large data sets, and so there is a need for an algorithm that efficiently produces consistently high-quality OTU assignments.Second Best (sounding): average neighbor
Notes: The alternative is to use our cluster.split command. 
These decreases in MCC value resulted in the formation of as many as 4.7 and 22.5% more OTUs, respectively, than were observed from the entire dataset. 
The use of the cluster splitting heuristic was probably not worth the loss in clustering quality. However, as datasets become larger, it may be necessary to use the heuristic to clustering the data into OTUs. 
Combine the clustering and count data into a human readable OTU table. The label is the percent similarity to required to make an OTU. What cutoff do you use? What level is a microbial genus? Species? Which do you want to look at?

Tutorial Command:
mothur > make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, label=0.03)
Notes: https://www.mothur.org/wiki/Make.shared


set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files)
make.shared(list=, count=, label=0.03?)
make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, label=0.03)
For Len298 Diff 1:
set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files)
make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, label=0.03)
Output File Names: 
/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff1.phylip.opti_mcc.shared
For Len298 Diff 2:
set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files)
make.shared(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, label=0.03)

set.dir(output=/home/micb405/Group12/Project3/clustering/shared_files/dgc)
make.shared(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, label=0.03)
Classify OTUs
Notes: You will notice that sequence AY457912 has a bootstrap value of 54% for the assignment to the Lachnospira pectinoschiza. This isn't much of a vote of confidence for this assignment. The current thinking seems to be to use a minimum cutoff of 60%. Mothur's default is set to a value of 80%:
SILVA vs. GreenGenes: https://www.researchgate.net/publication/315059693_SILVA_RDP_Greengenes_NCBI_and_OTT_-_how_do_these_taxonomies_compare
Silva: /home/micb405/data/project_3/databases/silva.nr_v128.align



Classify the sequences based on a database. Silva or GreenGenes? What bootstrap cutoff for confidence in your taxonomic assignment?

298diff2:
set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)
Doing: 400diff1
set.dir(output=/home/micb405/Group12/Project3/classifyseqs)
298diff1
set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)
400 diff1
set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)

400diff2
set.dir(output=/home/micb405/Group12/Project3/classifyseqs/cutoff60)
classify.seqs(fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.fasta, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, template=/home/micb405/data/project_3/databases/silva.nr_v128.align, taxonomy=/home/micb405/data/project_3/databases/silva.nr_v128.tax, cutoff=60)
Output File Names: 
/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy
/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.tax.summary


Use the sequence classification output to apply the resulting taxonomy to OTUs
298diff2 (80-80)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, cutoff=80, threshold=80, basis=otu) 
298diff2 (60-60)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff2.count_table, cutoff=60, threshold=60, basis=otu) 
298diff1 (80-80)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, cutoff=80, threshold=80, basis=otu) 
298diff1 (60-60)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8diff1.count_table, cutoff=60, threshold=60, basis=otu) 

400 diff 1 (80-80)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 1 (60-60)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 

400 diff 2 (80-80)
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, cutoff=80, threshold=80, basis=otu) 

Output File Names: 
/home/micb405/Group12/Project3/classifyotu/c80t80/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.0.03.cons.taxonomy
/home/micb405/Group12/Project3/classifyotu/c80t80/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.0.03.cons.tax.summary

400 diff 2 (60-60)

set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60)
classify.otu(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff2.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff2.count_table, cutoff=60, threshold=60, basis=otu) 

Output File Names: 
/home/micb405/Group12/Project3/classifyotu/c60t60/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.0.03.cons.taxonomy
/home/micb405/Group12/Project3/classifyotu/c60t60/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.0.03.cons.tax.summary


400 diff 1 (80-80)0.01
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80/0.01)
classify.otu(list=/home/micb405/Group12/Project3/clustering/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 1 (60-60)0.01
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60/0.01)
classify.otu(list=/home/micb405/Group12/Project3/clustering/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 
400 diff 1 (80-80)dgc
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c80t80/dgc)
classify.otu(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=80, threshold=80, basis=otu) 

400 diff 1 (60-60)dgc
set.dir(output=/home/micb405/Group12/Project3/classifyotu/c60t60/dgc)
classify.otu(list=/home/micb405/Group12/Project3/clustering/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.list, taxonomy=/home/micb405/Group12/Project3/classifyseqs/cutoff60/Saanich150m_t0.5k8_len400_diff1.nr_v128.wang.taxonomy, count=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.count_table, cutoff=60, threshold=60, basis=otu) 


Summarize OTU table
Summarize the data to get some descriptive statistics. There are many, many calc options found here.
summary.single(shared=, label=, calc=nseqs-sobs-coverage)
298diff1
set.dir(output=/home/micb405/Group12/Project3/summaryotu)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

298diff2
set.dir(output=/home/micb405/Group12/Project3/summaryotu)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8diff2.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)
400 diff 1
set.dir(output=/home/micb405/Group12/Project3/summaryotu)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)

400diff2
set.dir(output=/home/micb405/Group12/Project3/summaryotu)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)
Output File Names: 
/home/micb405/Group12/Project3/summaryotu/Saanich150m_t0.5k8_len400_diff2.phylip.opti_mcc.groups.summary
400 diff 1 0.01
set.dir(output=/home/micb405/Group12/Project3/summaryotu/0.01)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)
set.dir(output=/home/micb405/Group12/Project3/summaryotu/0.01)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/0.01/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared, calc=nseqs-sobs-coverage-shannon-chao)
400 diff 1 dgc
set.dir(output=/home/micb405/Group12/Project3/summaryotu/dgc)
summary.single(shared=/home/micb405/Group12/Project3/clustering/shared_files/dgc/Saanich150m_t0.5k8_len400_diff1.dgc.unique_list.shared, calc=nseqs-sobs-coverage-shannon-chao)










Tree:			
					
						
time mafft --maxiterate 1000 --auto /home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta >/home/micb405/Group12/Project3/mafft/400diff1.mfa 
					
				
			
		
time FastTree /home/micb405/Group12/Project3/mafft/400diff1.mfa 1>/home/micb405/Group12/Project3/mafft/tree 
					
				
Running the python Script:
Note: script is now located in /home/micb405/Group12/Project3/scripts/python_script/MICB405_project3/
Execute with:
 python script_name.py
OTU summary at genus level located in: /home/micb405/Group12/Project3/classifyotu/c60t60/ (for the c60t60)
0.01 located at: /home/micb405/Group12/Project3/classifyotu/c60t60/0.01
Dgc located at: /home/micb405/Group12/Project3/classifyotu/c60t60/dgc
There are 176 OTUs in your taxonomy file.
There are 81 distinct groups at the genus level in your taxonomy file.


C80t80 isnt working for me :(
Here are the files i’m using…. 
/home/micb405/Group12/Project3/classifyotu/c80t80/probs_f/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.0.03.cons.taxonomy
/home/micb405/Group12/Project3/clustering/shared_files/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.shared

Analysis Examples:
https://www.mothur.org/wiki/Agricultural_soil_community_analysis





FastTree -nt -gtr </home/micb405/Group12/Project3/mafft/400diff1.mfa 1>/home/micb405/Group12/Project3/mafft2/400_diff1.tree


set.dir(output=/home/micb405/Group12/Project3/mafft)
get.oturep(phylip=/home/micb405/Group12/Project3/otu_cluster/distance_matrix/len_400_diff_1/Saanich150m_t0.5k8_len400_diff1.phylip.dist, list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta)
time mafft --maxiterate 1000 --auto /home/micb405/Group12/Project3/mafft/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.0.03.rep.fasta >/home/micb405/Group12/Project3/mafft/400diff1tree.mfa 
					
				
time FastTree /home/micb405/Group12/Project3/mafft
/home/micb405/Group12/Project3/mafft/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.0.03.rep.fasta 1>/home/micb405/Group12/Project3/mafft/kk.tree

bin.seqs(list=/home/micb405/Group12/Project3/clustering/Saanich150m_t0.5k8_len400_diff1.phylip.opti_mcc.list, fasta=/home/micb405/Group12/Project3/finalfasta/Saanich150m_t0.5k8_len400_diff1.fasta )





