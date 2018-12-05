##Super-Scaffold_84 Validation

We are looking at a subsection of potential cross over sites. I picked Super-Scaffold_84 because it has a weird 
inverted repeat structure that seems to be created by the joining of syntenic homeologous regions to a middle
syntentic block. 

![Super-Scaffold_84 Structure]



This is potentially really cool. One question I have before we get to excited is how much can we trust this assembly.
Previous work with mate-pair reads suggest that scaffold of all cross-over(xo) regions is trusty. In that mate-pair 
reads mapped with bowtie2 --local provide continous links across a xo regions investigated. MP mapping is not insanely
deep, compared to say the mate-pairs you would pull down on plastid. I wanted to get a finer grain view of some of these
junctions. Particularly, 84 -> 153 and 153->84  syntenic blocks. 

To this end I used `bowtie2` to map across Super-Scaffold_84.

```bash 
 nohup bowtie2 -x ss84.fa \
 --threads 15 \
 --no-unal \
 --very-fast \
 --trim5 10 \
 --trim3 10 \
 -1 /home/biolinux/scratch/Eleusine_Synteny_Proj/data/raw_wgs_reads/DRR095893_1.fastq \
 -2 /home/biolinux/scratch/Eleusine_Synteny_Proj/data/raw_wgs_reads/DRR095893_2.fastq \
 -S ss84_DRR095893_vf53trim.sam &>map.DRR095893.nohup & 
 
ndh0004@beta:/home/biolinux/scratch/Eleusine_Synteny_Proj/data/primer_dev/ss84$ bowtie2 --version
/automnt/opt/bin/../bowtie2/default/bowtie2-align-s version 2.2.9
64-bit
Built on localhost.localdomain
Thu Apr 21 18:36:37 EDT 2016
Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-54)
Options: -O3 -m64 -msse2  -funroll-loops -g3 -DPOPCNT_CAPABILITY
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}


samtools view -hb ss84_DRR095893_vf53trim.sam | samtools sort - ss84_DRR095893_vf53trim 
samtools index ss84_DRR095893_vf53trim.bam


```
 gff loaded from bog_sept28*gff3.
 
I examined junctions in tablet. I found that there was tight
paired end read mapping across each of these regions.
![Example Junction][ss84_junc]

  For example when we look at a blast comparison of this junction we see differences.
  
![Blast Result][bres]


This is encouraging because we found by blasting around the ends of this junction that there appear to repeated
sequence surrounding this junction. Though not perfectly repeated. 
 
 ```bash 
 example blasts
 # BLASTN 2.2.31+
# Query: Super-Scaffold_158:507103-509103
# Database: ../Ecor_PR202_scaffolds.fa
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 7 hits found
Super-Scaffold_158:507103-509103	Super-Scaffold_158	100.00	2001	0	0	1	2001	507103	509103	0.0	3696
Super-Scaffold_158:507103-509103	Super-Scaffold_84	93.98	1013	49	7	461	1469	3413571	3414575	0.0	1522
Super-Scaffold_158:507103-509103	Super-Scaffold_84	89.76	332	28	5	94	422	3411353	3411681	7e-115	420
Super-Scaffold_158:507103-509103	Super-Scaffold_84	85.38	171	3	7	1843	2001	3414909	3415069	6e-36	158
Super-Scaffold_158:507103-509103	Super-Scaffold_84	92.86	98	7	0	1	98	3411067	3411164	2e-31	143
Super-Scaffold_158:507103-509103	Super-Scaffold_84	88.37	86	6	2	1497	1579	3414741	3414825	1e-18	100
Super-Scaffold_158:507103-509103	Super-Scaffold_16	93.75	48	1	1	1796	1843	2504462	2504417	8e-10	71.3
# BLAST processed 1 queries
 
 
 ```
 
 
 

 
[Super-Scaffold_84 Structure]:ss84_structure.png
[ss84_junc]:tablet_ss84_junc.png
[bres]:ss84_ss158_junc.png