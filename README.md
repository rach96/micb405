# micb405: Bash Scripts/Commands
Project 2: Metagenomics

**Task:**

Reference the bash script(s) containing all commands. These are to be either hosted on a group member's GitHub page and/or emailed to Connor and I will upload them to the MICB405-Metagenomics GitHub page under student_scripts/Group*/.

**Workflow Diagram:**

![screen shot 2017-11-21 at 6 47 59 pm](https://user-images.githubusercontent.com/25336570/33107421-0d96b142-ceec-11e7-8570-591d05fec779.png)

**Directory Structure:**

![screen shot 2017-11-21 at 7 18 57 pm](https://user-images.githubusercontent.com/25336570/33108273-48a16490-cef0-11e7-8bcb-82dfaa1f09a6.png)

**FastQC Commands:** 
```
fastqc --threads 2 -o /home/micb405/Group12/Project2/FastQC_Output/ \
/home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz 

fastqc --threads 2 -o /home/micb405/Group12/Project2/FastQC_Output/ \
/home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz
```

**MEGAHIT (Li et al. 2015) Commands:**

```
nohup megahit -1 /home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz -2 \
/home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz --k-min 27 --k-max 147 \
--k-step 20 --min-contig-len 1000 -m 0.07 -t 2 --out-dir \
/home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m & 
```

**MaxBin (Version 2.0; Wu et al. 2014) Commands:**

```
export PATH=/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:/usr/local/sbin: \
/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/micb405/resources/project_2/FragGeneScan1.30
```
```
nohup perl5.26.0 /home/micb405/resources/project_2/MaxBin-2.2.4/run_MaxBin.pl -contig \ 
/home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa -reads \
/home/micb405/data/project_2/SI072_LV_150m_DNA_R1.fastq.gz -reads2 \
/home/micb405/data/project_2/SI072_LV_150m_DNA_R2.fastq.gz -out myout -thread 2 -plotmarker &
```

**CheckM (Parks et al. 2014) Commands:**

```
checkm lineage_wf --tab_table -x .fasta --threads 4 --pplacer_threads 4 /home/micb405/Group12/Project2/MaxBin_output/ \
/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkm_output/ \  
>/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkM_stdout.tsv
```
* Note: This command was run by Connor

**Exporting the CheckM Data (From Connor's Workflow):**

```
awk -F"\t" '{ if ($12>10 && $13<5) print $0 }' /home/micb405/Group12/Project2/checkM_output/Group12_checkM_stdout_file.tsv > \ 
/home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv
```

**Taxonomic Classification Using Mash: (From Connor's Workflow)**

**1) BASH Script Using the RefSeq (Pruitt et al. 2007) Database**

```
#!/bin/bash                                                                                                                      
while read line
do
bin=$( echo $line | awk '{ print $1 }')
sid=$( echo $bin | awk -F. '{ print $1 }')
if [ -f /home/micb405/Group12/Project2/MaxBin_Good/myout.006.fasta ]
    then
    mash dist -v 1E-8 /home/micb405/resources/project_2/refseq.genomes.k21s1000.msh \
    /home/micb405/Group12/Project2/MaxBin_Good/myout.001.fasta
fi
done</home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv \
>/home/micb405/Group12/Project2/tables/RefSeq_Mash_output_001.tsv
```
* Note: This command was repeated for bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)

**2) After the output .tsv files were generated, the highest quality annotations were extracted using the command (from Connor's Workflow):**

```
cat RefSeq_Mash_output.tsv Saanich_Mash_output.tsv | sort -t$'\t' -k2,2 | \
awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' >Mash_classifications.BEST.tsv
```

**Taxonomic Classification Using LAST and the Silva (Quast et al. 2013) Database:**

**1) Command Using the Silva Database**

```
while read line; do bin=$( echo $line | awk '{ print $1 }'); sid=$( echo $bin | awk -F. \
'{ print $1 }'); if [ -f /home/micb405/Group12/Project2/MaxBin_Good/myout.001.fasta ]; \
then best_hit=$(lastal -f TAB -P 4 /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva \ /home/micb405/Group12/Project2/MaxBin_Good/myout.001.fasta | grep -v "^#" | head -1); \
echo $bin,$sid,$best_hit | sed 's/,\| /\t/g'; fi; \
done</home/micb405/Group12/Project2/tables/GT10Complete_LT5Contam_MAGs_checkM.tsv \ >/home/micb405/Group12/Project2/LAST_tables/LAST_SILVA_alignments_001.BEST.tsv
```
**2) Editing the output file to display taxonomic information**

```
while read line; do accession=$( echo $line | awk '{ print $4 }'); bin=$( echo $line | awk \
'{ print $1 }' ); if [ ! -z $accession ]; then last_hit=$( grep "$accession" \ /home/micb405/resources/project_2/SILVA_128_SSURef_taxa_headers.txt | awk \ 
'{ $1=""; print $0 }'); echo $bin,$last_hit; fi; \ 
done</home/micb405/Group12/Project2/LAST_tables/LAST_SILVA_alignments_069.BEST.tsv >/home/micb405/Group12/Pr
```
* Note: This command was repeated for bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%)


**Combining the MASH and LAST Data: Generating the TSV File**

```
cat RefSeq_Mash_output_001.tsv RefSeq_Mash_output_006.tsv RefSeq_Mash_output_007.tsv \
RefSeq_Mash_output_019.tsv RefSeq_Mash_output_021.tsv RefSeq_Mash_output_028.tsv \
RefSeq_Mash_output_046.tsv RefSeq_Mash_output_058.tsv RefSeq_Mash_output_065.tsv \
RefSeq_Mash_output_069.tsv Saanich_Mash_output_001.tsv Saanich_Mash_output_006.tsv \ 
Saanich_Mash_output_007.tsv Saanich_Mash_output_019.tsv Saanich_Mash_output_021.tsv \
Saanich_Mash_output_028.tsv Saanich_Mash_output_046.tsv Saanich_Mash_output_058.tsv \
Saanich_Mash_output_065.tsv Saanich_Mash_output_069.tsv | sort -t$'\t' -k2,2 | \
awk '{ if(!x[$2]++) {print $0; dist=($3-1)} else { if($3<dist) print $0} }' > \ 
/home/micb405/Group12/Project2/tables/Mash_classifications.BEST.tsv
```

**PROKKA (Prokka, 2014) Commands:**

```
prokka --prefix myout.001 /home/micb405/Group12/Project2/MaxBin_output/myout.001.fasta 
```
* Note: This command was repeated for bins: 6,7,9,19,21,24,28,46,58,65,68,69 (met the threshold of completeness > 10% and contamination < 5%) 

**RPKM (Mortazavi et al. 2008) Commands:**

```
nohup bwa index /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa &
```

```
nohup bwa mem -t 4 /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa \
/home/dtruong/SI072_LV_150m_DNA_R1.fastq.gz /home/dtruong/SI072_LV_150m_DNA_R2.fastq.gz \
1>/home/micb405/Group12/Project2/BWA_output/SI072_LV_150m_DNA.sam \
2>/home/micb405/Group12/Project2/BWA_output/SI072_LV_150m_DNA.bwa.stderr &
```

```
/home/micb405/resources/project_2/rpkm -c /home/micb405/Group12/Project2/MEGAHIT/SI072_LV_150m/final.contigs.fa \
-a /home/micb405/Group12/Project2/BWA_output/SI072_LV_150m_DNA.sam -o \
/home/micb405/Group12/Project2/RPKM/SI072_LV_150m_DNA_RPKM.csv
```

```
ls /home/micb405/Group12/Project2/MaxBin_Good/*fasta >mag_list.txt
```

```
/home/micb405/resources/project_2/find_mag_rpkm_average.py -l \
/home/micb405/Group12/Project2/MaxBin_Good/mag_list.txt -r \
/home/micb405/Group12/Project2/RPKM/SI072_LV_150m_DNA_RPKM.csv \
-o /home/micb405/Group12/Project2/RPKM/SI072_LV_150m_MAG_RPKM.csv
```

## References:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Glöckner FO, Yilmaz P, Quast C, Gerken J, Beccati A, Ciuprina A, Bruns G, Yarza P, Peplies J, Westram R, Ludwig W (2017) 25 years of serving the community with ribosomal RNA gene reference databases and tools. J. Biotechnol. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, doi: 10.1093/bioinformatics/btv033 [PMID: 25609793]. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.
    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Mortazavi, A., Williams, B.A., McCue, K., Schaeffer, L., Wold, B. (2008). Mapping and quantifying mammalian transcriptomes by RNA-Seq. Nature Methods 5, 621-628 doi:10.1038/nmeth.1226
    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol. 2016 Jun 20;17(1):132. doi: 10.1186/s13059-016-0997-x.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2014. Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043-1055.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pruitt, K. D., Tatusova, T., & Maglott, D. R. (2007). NCBI reference sequences (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins. Nucleic Acids Research, 35(Database issue), D61–D65. http://doi.org/10.1093/nar/gkl842

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Seemann T. Prokka: rapid prokaryotic genome annotation Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Wu YW, Tang YH, Tringe SG, Simmons BA, and Singer SW, "MaxBin: an automated binning method to recover individual genomes from metagenomes using an expectation-maximization algorithm", Microbiome, 2:26, 2014.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Wu YW, Simmons BA, and Singer SW, "MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets", Bioinformatics, 32(4): 605-607, 2016.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648


Still Don't Have:
* FASTQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* LAST: http://last.cbrc.jp/doc/last-papers.html
