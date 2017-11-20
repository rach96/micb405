# micb405: Bash Scripts/Commands
Project 2: Metagenomics

**Task:**

Reference the bash script(s) containing all commands. These are to be either hosted on a group member's GitHub page and/or emailed to Connor and I will upload them to the MICB405-Metagenomics GitHub page under student_scripts/Group*/.

**FastQC Commands:** 

SI072_LV_150m_DNA_R1.fastq.gz
SI072_LV_150m_DNA_R2.fastq.gz
fastqc --threads 2 -o /home/micb405/Group12/Project2/FastQC_Output/ /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz
fastqc --threads 2 -o /home/micb405/Group12/Project2/FastQC_Output/ /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz

**MEGAHIT Commands:**

nohup megahit -1 /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz -2 /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz --k-min 27 --k-max 147 --k-step 20 --min-contig-len 1000 -m 0.07 -t 2 --out-dir /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m &

**MaxBin2 Commands:**

export PATH=/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/micb405/resources/project_2/FragGeneScan1.30

nohup perl5.26.0 /home/micb405/resources/project_2/MaxBin-2.2.4/run_MaxBin.pl -contig /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa -reads /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz -reads2 /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz -out myout -thread 2 -plotmarker &

**CheckM Commands:**

checkm lineage_wf --tab_table -x .fasta --threads 4 --pplacer_threads 4 $BIN_DIR \
/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkm_output/ >/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkM_stdout.tsv
* Note: This command was run by Connor

**Exporting the CheckM Data:**

awk -F"\t" '{ if ($12>10 && $13<5) print $0 }' /home/micb405/Group12/Project2/checkM_output/Group12_checkM_stdout_file.tsv > /home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv

**MASH Commands:**

mash dist /home/micb405/resources/project_2/refseq.genomes.k21s1000.msh /home/micb405/Group12/Project2/MaxBin_output/myout.001.fasta > /home/micb405/Group12/Project2/Mash_output/myout.001.mash
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

**Exporting the MASH Data: BASH Scripts**

**1) BASH Script Using the RefSeq Database**

#!/bin/bash                                                                                                                      
while read line
do
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f /home/micb405/Group12/Project2/MaxBin_Good/myout.006.fasta ]
    then
    mash dist -v 1E-8 /home/micb405/resources/project_2/refseq.genomes.k21s1000.msh /home/micb405/Group12/Project2/MaxBin_Good/m\
yout.001.fasta
fi
done</home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv >/home/micb405/Group12/Project2/tables/RefSeq_\
Mash_output_001.tsv
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

**2) BASH Script Using the 16S rRNA Database**

#!/bin/bash                                                                                                                      

while read line
do
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f /home/micb405/Group12/Project2/MaxBin_Good/myout.006.fasta ]
    then
    mash dist -v 1E-8 /home/micb405/resources/project_2/Saanich_QCd_SAGs_k21s1000.sig.msh /home/micb405/Group12/Project2/MaxBin_\
Good/myout.001.fasta
fi
done</home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv >/home/micb405/Group12/Project2/tables/Saanich\
_Mash_output_001.tsv
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

**Exporting the MASH Data: Generating the TSV File**

cat RefSeq_Mash_output.tsv Saanich_Mash_output.tsv | sort -t$'\t' -k2,2 | awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' >Mash_classifications.BEST.tsv

cat RefSeq_Mash_output_001.tsv RefSeq_Mash_output_006.tsv RefSeq_Mash_output_007.tsv RefSeq_Mash_output_019.tsv RefSeq_Mash_output_021.tsv RefSeq_Mash_output_028.tsv
RefSeq_Mash_output_046.tsv RefSeq_Mash_output_058.tsv RefSeq_Mash_output_065.tsv RefSeq_Mash_output_069.tsv Saanich_Mash_output_001.tsv Saanich_Mash_output_006.tsv Saanich_Mash_output_007.tsv Saanich_Mash_output_019.tsv Saanich_Mash_output_021.tsv Saanich_Mash_output_028.tsv Saanich_Mash_output_046.tsv Saanich_Mash_output_058.tsv Saanich_Mash_output_065.tsv Saanich_Mash_output_069.tsv | sort -t$'\t' -k2,2 | awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' > /home/micb405/Group12/Project2/tables/Mash_classifications.BEST.tsv

**LAST Commands: Running lastal**

lastal -f TAB /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva /home/micb405/Group12/Project2/MaxBin_output/myout.001.fasta >/home/micb405/Group12/Project2/LAST_output/myout.001.tab
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%) 

**Exporting the LAST Data: Generating the TSV File**

while read line; do bin=$( echo $line | awk '{ print $1 }'); sid=$( echo $bin | awk -F. '{ print $1 }'); if [ -f /home/micb405/Group12/Project2/MaxBin_Good/myout.001.fasta ]; then best_hit=$(lastal -f TAB -P 4 /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva /home/micb405/Group12/Project2/MaxBin_Good/myout.069.fasta | grep -v "^#" | head -1); echo $bin,$sid,$best_hit | sed 's/,\| /\t/g'; fi; done</home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv >/home/micb405/Group12/Project2/LAST_tables/LAST_SILVA_alignments_001.BEST.tsv
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

while read line; do accession=$( echo $line | awk '{ print $4 }'); bin=$( echo $line | awk '{ print $1 }' ); if [ ! -z $accession ]; then last_hit=$( grep "$accession" /home/micb405/resources/project_2/SILVA_128_SSURef_taxa_headers.txt | awk '{ $1=""; print $0 }'); echo $bin,$last_hit; fi; done</home/micb405/Group12/Project2/LAST_tables/LAST_SILVA_alignments_001.BEST.tsv >/home/micb405/Group12/Project2/LAST_tables/LAST_SILVA_classifications_001.BEST.csv
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

**PROKKA Commands:**

prokka --prefix myout.001 /home/micb405/Group12/Project2/MaxBin_output/myout.001.fasta 
* Note: This command was repeated for Bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%) 

**RPKM Commands:**

nohup bwa index /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa &

Copied SI072_LV_150m_DNA_R1.fastq.gz and SI072_LV_150m_DNA_R2.fastq.gz into /home/dtruong/ directory and re-ran the command as

nohup bwa mem -t 4 /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa
/home/dtruong/SI072_LV_150m_DNA_R1.fastq.gz /home/dtruong/SI072_LV_150m_DNA_R2.fastq.gz
1>/home/micb405/Group12/Project2/BWA_output/SI072_LV_150m_DNA.sam 2>/home/micb405/Group12/Project2/BWA_output/SI072_LV_150m_DNA.bwa.stderr &

/home/micb405/resources/project_2/rpkm -c /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa
-a /home/micb405/Group12/Project2/BWA_output/SI072_LV_150m_DNA.sam -o /home/micb405/Group12/Project2/RPKM/SI072_LV_150m_DNA_RPKM.csv

ls /home/micb405/Group12/Project2/MaxBin_Good/*fasta >mag_list.txt

/home/micb405/resources/project_2/find_mag_rpkm_average.py -l /home/micb405/Group12/Project2/MaxBin_Good/mag_list.txt -r /home/micb405/Group12/Project2/RPKM/SI072_LV_150m_DNA_RPKM.csv
-o /home/micb405/Group12/Project2/RPKM/SI072_LV_150m_MAG_RPKM.csv




