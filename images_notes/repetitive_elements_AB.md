## Repeat calling

1 started with repeat modeler
```bash 
: cat /scratch/ndh0004/databases/Ecor_genome/*sh
#!/bin/bash 
# ----------------QSUB Parameters----------------- #
##choose queue
####PBS -q
##list - node are nodes: ppn are cpus per node: walltime=walltime
#PBS -l nodes=2:ppn=5,mem=204gb,walltime=40:00:00
##email
#PBS -M <your-email_here>
##send email abort; begin; end
#PBS -m ae
##job name
#PBS -N bt2_ind
##combine standard out and standard error
#PBS -j oe
##queue reservation: this may change keep an eye on it!

# ----------------Load Modules-------------------- #

cd /scratch/ndh0004/databases/Ecor_genome

module load repeatmodeler

RepeatModeler -engine ncbi -database Ecor_PR202_scaffolds -pa 10

```

2 Ran modeled repeats against each A and B regions B script shown flags are identical to A

```cat repeatMasker.sh 
#!/bin/bash 
# ----------------QSUB Parameters----------------- #
##choose queue
####PBS -q
##list - node are nodes: ppn are cpus per node: walltime=walltime
#PBS -l nodes=1:ppn=1,walltime=70:00:00
##email
##PBS -M <your-email_here>
##send email abort; begin; end
#PBS -m ae
##job name
#PBS -N bt2_ind
##combine standard out and standard error
#PBS -j oe
##queue reservation: this may change keep an eye on it!

# ----------------Load Modules-------------------- #

cd /scratch/ndh0004/repetitive_element_coracana/B

module load repeatmasker

RepeatMasker -s -lib consensi.fa.classified -e ncbi B.fasta

[ndh0004@hopper-login](B)[14:26]: pwd 
/scratch/ndh0004/repetitive_element_coracana/B
```

3 get [AB].fasta.tbl
```bash
: cat */*tbl
==================================================
file name: A.fasta                  
sequences:            68
total length:  128821903 bp  (126120816 bp excl N/X-runs)
GC level:         44.05 %
bases masked:   53704138 bp ( 41.69 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:             1374       276806 bp    0.21 %
      ALUs            0            0 bp    0.00 %
      MIRs            0            0 bp    0.00 %

LINEs:             4911      2535727 bp    1.97 %
      LINE1        4136      2210551 bp    1.72 %
      LINE2          18         9976 bp    0.01 %
      L3/CR1          0            0 bp    0.00 %

LTR elements:     24930     37237186 bp   28.91 %
      ERVL            0            0 bp    0.00 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI    110        15076 bp    0.01 %
      ERV_classII     0            0 bp    0.00 %

DNA elements:     22920      7682996 bp    5.96 %
     hAT-Charlie     50         6362 bp    0.00 %
     TcMar-Tigger     0            0 bp    0.00 %

Unclassified:     13285      4864215 bp    3.78 %

Total interspersed repeats: 52596930 bp   40.83 %


Small RNA:         1479       332446 bp    0.26 %

Satellites:           0            0 bp    0.00 %
Simple repeats:   21350       933800 bp    0.72 %
Low complexity:    2395       119456 bp    0.09 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element
                                                      

The query species was assumed to be homo          
RepeatMasker version open-4.0.6 , sensitive mode
                                 
run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829
==================================================
file name: B.fasta                  
sequences:            82
total length:  165069540 bp  (154448130 bp excl N/X-runs)
GC level:         43.68 %
bases masked:   75739598 bp ( 45.88 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:             1772       322755 bp    0.20 %
      ALUs            0            0 bp    0.00 %
      MIRs            0            0 bp    0.00 %

LINEs:             5416      2903384 bp    1.76 %
      LINE1        4282      2319530 bp    1.41 %
      LINE2          51        19778 bp    0.01 %
      L3/CR1          0            0 bp    0.00 %

LTR elements:     33307     56537939 bp   34.25 %
      ERVL            0            0 bp    0.00 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI    141        21333 bp    0.01 %
      ERV_classII     0            0 bp    0.00 %

DNA elements:     26746      9314705 bp    5.64 %
     hAT-Charlie     54         6828 bp    0.00 %
     TcMar-Tigger     0            0 bp    0.00 %

Unclassified:     15752      5488288 bp    3.32 %

Total interspersed repeats: 74567071 bp   45.17 %


Small RNA:         1853       352446 bp    0.21 %

Satellites:           0            0 bp    0.00 %
Simple repeats:   22780      1019634 bp    0.62 %
Low complexity:    2772       136175 bp    0.08 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element
                                                      

The query species was assumed to be homo          
RepeatMasker version open-4.0.6 , sensitive mode
                                 
run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829

```

4 create a set of bed files to look at coverage

```bash 
awk ' OFS="\t" {gsub("/","_", $11);  print $5,$6,$7 > $11".bed" }' B.fasta.out
```

5 create uniform sliding windows

```bash 
# start with faidx create a genome length list 
awk '{print $1"#"$2}' B.fasta.fai > genome.lengths
# then create silding windows

awk ' OFS="\t" {gsub("/","_", $11);  print $5,$6,$7 > $11".bed" }' B.fasta.out

# not to bed files with coverage

for X in $( cat genome.lengths ); do  len=$( echo $X | awk  -F'#' '{print $2}' ) ; genome=$( echo $X |  awk  -F'#' '{print $1}' ) ; start=0 ; while [ $(( $start + 100000 )) -lt $len ] ; do printf "${genome}\t${start}\t$(( $start + 100000))\n" ; start=$(( $start + 20000 )) ; done   ; done  > B_sliding_window.bed

# now get the coverage
for X in   $( ls *bed | egrep -v sliding ) ;  do bedtools  coverage  -a B_sliding_window.bed -b $X > coverage_${X} ;echo $X ;  done

# Now let's put the all into 1 summary file

awk -F"\t" ' BEGIN {print"contig start stop"}; {print $1,$2,$3 }' coverage_DNA.bed > sliding_window_coverarage.csv && for X in coverage_*bed ; do awk -F"\t" -vT=$( echo $X | sed 's/.bed//;  s/coverage_//g' )  ' BEGIN {print"count_"T" bp_"T}; {print $4,$5 }' $X | paste sliding_window_coverarage.csv - | sed 's/\t/ /g'>tmp ; mv tmp sliding_window_coverarage.csv ;done

mv sliding_window_coverarage.csv  B_sliding_window_coverarage.csv 
# this was done with both A and B
```


## Genome Painting with A vs B Type I retrotransposons


1 start by getting families only represented in both genomes. We are assuming that for this case that a majoritiy of the
time there is no to little movement between A and B genomes. 

```bash

# let's take sections longer than 200bp

: awk '($11 ~ /LTR/ || $11 ~ /SINE/ || $11 ~ /SINE/) && $7-$6 >200 {print $10 }' ../B/B.fasta.out | sort |uniq > B_uniq_subfamily.list 

: awk '($11 ~ /LTR/ || $11 ~ /SINE/ || $11 ~ /SINE/) && $7-$6 >200 {print $10 }' ../A/A.fasta.out | sort |uniq > A_uniq_subfamily.list 

#

[13:59]: awk '$10 != "rnd-4_family-2017" && $10 != "rnd-5_family-2034" && $10 != "rnd-5_family-80" && $1 != "rnd-6_family-3835" && $7-$6 >200  {OFS="\t"; print $5,$6,$7 }' ../A/A.fasta.out >Aretrotransposon.bed

[14:01]: awk '$10 != "rnd-4_family-2017" && $10 != "rnd-5_family-2034" && $10 != "rnd-5_family-80" && $1 != "rnd-6_family-3835" && $7-$6 >200  {OFS="\t"; print $5,$6,$7 }' ../B/B.fasta.out >Bretrotransposon.bed
 qsub -I 
 
 [ndh0004@node179 for_genome_painting]$ bedtools getfasta -fi ../A/A.fasta -bed Aretrotransposon.bed | sed 's/>/>A_/' >retro.fasta
[ndh0004@node179 for_genome_painting]$ bedtools getfasta -fi ../B/B.fasta -bed Bretrotransposon.bed | sed 's/>/>B_/' >>retro.fasta


[ndh0004@node179 for_genome_painting]$ grep -c \> retro.fasta 
105091

[ndh0004@node179 for_genome_painting]$ module load bowtie2
[ndh0004@node179 for_genome_painting]$ bowtie2-build -q retro.fasta retro.fasta 


```

1 now lets get A and B retrotransposable reads


```bash

#!/bin/bash 
# ----------------QSUB Parameters----------------- #
##choose queue
####PBS -q
##list - node are nodes: ppn are cpus per node: walltime=walltime
#PBS -l nodes=5:ppn=4,walltime=70:00:00
##email
##PBS -M <your-email_here>
##send email abort; begin; end
#PBS -m ae
##job name
#PBS -N bt2_ind
##combine standard out and standard error
#PBS -j oe
##queue reservation: this may change keep an eye on it!

# ----------------Load Modules-------------------- #

module load bowtie2
module load samtools 
cd /scratch/ndh0004/repetitive_element_coracana/for_genome_painting
fq_L=/scratch/ndh0004/repetitive_element_coracana/for_genome_painting/DRR095893_1.fastq  
fq_R=/scratch/ndh0004/repetitive_element_coracana/for_genome_painting/DRR095893_2.fastq
ref=/scratch/ndh0004/repetitive_element_coracana/for_genome_painting/retro.fasta

threads=20
bowtie2 -x ${ref}  \
                --local \
                --threads $threads \
                -1 $fq_L \
                -2 $fq_R  | \
                samtools view -h -b - | samtools sort -n - -o tmp.bam

A_contigs=$( samtools idxstas tmp.bam | egrep '^A_' | awk '$1 != "*" {print $1}' )
B_contigs=$( samtools idxstas tmp.bam | egrep '^B_' | awk '$1 != "*" {print $1}' )

samtools index tmp.bam
samtools view -@${threads} -h -F 8 -f 4 tmp.bam > tmp_mapped.sam
samtools view -@${threads}  -f 8 -F 4 tmp.bam >> tmp_mapped.sam
samtools view -@5${threads}  -F 8 -F 4 tmp.bam >> tmp_mapped.sam
samtools view -hb tmp_mapped.sam | samtools sort -@${threads} -n - -o mapped.bam




```