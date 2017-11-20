# micb405: Bash Scripts/Commands
Project 2: Metagenomics

Task:

Reference the bash script(s) containing all commands. These are to be either hosted on a group member's GitHub page and/or emailed to Connor and I will upload them to the MICB405-Metagenomics GitHub page under student_scripts/Group*/.

FastQC Commands: 

SI072_LV_150m_DNA_R1.fastq.gz
SI072_LV_150m_DNA_R2.fastq.gz
fastqc --threads 2 -o /home/micb405/Group12/Project2/FastQC_Output/ /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz
fastqc --threads 2 -o /home/micb405/Group12/Project2/FastQC_Output/ /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz

MEGAHIT Commands:

nohup megahit -1 /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz -2 /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz --k-min 27 --k-max 147 --k-step 20 --min-contig-len 1000 -m 0.07 -t 2 --out-dir /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m &

MaxBin2 Commands:

export PATH=/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/micb405/resources/project_2/FragGeneScan1.30

nohup perl5.26.0 /home/micb405/resources/project_2/MaxBin-2.2.4/run_MaxBin.pl -contig /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa -reads /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz -reads2 /home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz -out myout -thread 2 -plotmarker &

CheckM Commands:

checkm lineage_wf --tab_table -x .fasta --threads 4 --pplacer_threads 4 $BIN_DIR \
/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkm_output/ >/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkM_stdout.tsv
* Note: This command was run by Connor

MASH Commands:

mash dist /home/micb405/resources/project_2/refseq.genomes.k21s1000.msh /home/micb405/Group12/Project2/MaxBin_output/myout.058.fasta > /home/micb405/Group12/Project2/Mash_output/myout.058.mash
* Note: This command was repeated for Bins: 1,6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

LAST Commands:






