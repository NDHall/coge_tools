###Final calls for xo regions.
What we want to do first is to validate the assemblies in the regions of interest.
This can be useful for Hatekayama et al. (2018) data and for Hittalmani et al. (2017).
 
First we want to begin by getting a list of all AB contigs. 
I have produced notes on the underlying calls we have made.

![notes_on_AB](:https://github.com/NDHall/coge_tools/blob/master/images_notes/process_for_calling_AvsB.md)
![notes_on_dc_ksbagofgenes](:https://github.com/NDHall/coge_tools/blob/master/images_notes/notes_dc_ks_bagofgenes.md)
![notes_on_xo_sep24](:https://github.com/NDHall/coge_tools/blob/master/images_notes/notes_on_xo_sep24.md)



```bash
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep '^AB' bog_sep28_qaccalled_chrom_calls_all.tsv | sort -k2 -k3
AB	provisional	12	1	Super-Scaffold_70	104	115	26	20
AB	provisional	12	2	Super-Scaffold_95	187	24	20	12
AB	provisional	13	1	Super-Scaffold_10	114	126	0	1
AB	provisional	14	3	Super-Scaffold_50	220	65	9	15
AB	provisional	20	1	Super-Scaffold_39	65	292	12	19
AB	provisional	20	4	Super-Scaffold_84	225	250	30	37
AB	provisional	2	1	Super-Scaffold_2765	36	18	5	6
AB	provisional	3	1	scaffold306_size1284463_obj	38	44	5	5
AB	provisional	3	2	Super-Scaffold_202	5	17	9	2
AB	provisional	4	1	Super-Scaffold_147	12	50	11	4
AB	provisional	5	1	Super-Scaffold_210	65	13	4	11
AB	provisional	5	1	Super-Scaffold_76	37	13	1	0
AB	provisional	5	2	Super-Scaffold_291	29	34	7	3
AB	provisional	6	1	Super-Scaffold_266	24	71	1	5
AB	provisional	6	2	Super-Scaffold_158	71	49	26	19
AB	provisional	7	1	Super-Scaffold_25	126	27	4	1
AB	strong	10	0	Super-Scaffold_253	150	20	0	0
AB	strong	12	0	Super-Scaffold_128	27	297	0	0
AB	strong	12	0	Super-Scaffold_2412	35	225	0	0
AB	strong	14	0	Super-Scaffold_12	65	220	0	0
AB	strong	2	0	scaffold591_size730508_subseq_1:450970_obj	11	26	0	0
AB	strong	2	0	scaffold93_size1230182_subseq_1:987274_obj	14	15	0	0
AB	strong	4	0	Super-Scaffold_105	21	7	0	0
AB	strong	4	0	Super-Scaffold_1173	43	41	0	0
AB	strong	4	0	Super-Scaffold_13	41	43	0	0
AB	strong	6	0	Super-Scaffold_146	19	72	0	0
AB	strong	7	0	Super-Scaffold_307	25	83	0	0
AB	strong	9	0	Super-Scaffold_311	134	31	0	0

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep '^AB' bog_sep28_qaccalled_chrom_calls_all.tsv | sort -k2 -k3 | wc -l 
28
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ 
```

28 entries scaffolds to deal with. Hatekeyama et al. Have mp libraries 
with 20 kb inserts. We can now take the output of this code as a list.
get bed regions and determine where the cross overs take place. 

```bash
# take chrom name
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep '^AB' bog_sep28_qaccalled_chrom_calls_all.tsv | \
awk '{print $5}' | sort    > ab_scaf_list.txt
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f ab_scaf_list.txt  bog_sep28_qaccalled_region_calls.tsv | awk '{print $2}'| sort | uniq -c | wc -l 
28
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f ab_scaf_list.txt  bog_sep28_qaccalled_region_calls.tsv | awk '{print $2}'| sort | uniq -c | sort -n 
      2 scaffold591_size730508_subseq_1:450970_obj
      2 scaffold93_size1230182_subseq_1:987274_obj
      3 Super-Scaffold_2765
      4 scaffold306_size1284463_obj
      4 Super-Scaffold_105
      4 Super-Scaffold_1173
      4 Super-Scaffold_13
      5 Super-Scaffold_147
      5 Super-Scaffold_202
      6 Super-Scaffold_146
      6 Super-Scaffold_210
      6 Super-Scaffold_76
      7 Super-Scaffold_266
      7 Super-Scaffold_291
      7 Super-Scaffold_307
      8 Super-Scaffold_158
      8 Super-Scaffold_25
      9 Super-Scaffold_311
     10 Super-Scaffold_253
     12 Super-Scaffold_128
     12 Super-Scaffold_2412
     13 Super-Scaffold_70
     14 Super-Scaffold_10
     14 Super-Scaffold_12
     14 Super-Scaffold_95
     17 Super-Scaffold_50
     21 Super-Scaffold_39
     24 Super-Scaffold_84
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f ab_scaf_list.txt  bog_sep28_qaccalled_region_calls.tsv | awk ' OFS="\t" {print $2,$3,$4,$1}'| sort  |head -n5
scaffold306_size1284463_obj	226643	788390	B
scaffold306_size1284463_obj	724	84454	B
scaffold306_size1284463_obj	785065	900798	ambig
scaffold306_size1284463_obj	964484	1279227	A
scaffold591_size730508_subseq_1:450970_obj	260001	450970	B
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ ndh0004@IorekByrnison:~/code/coge_tools/data_out$ bedtools sort -i xo_chroms_oct1_2018.bed >sorted_xo_chroms_oct1_2018.bed


```
 Used `python /home/ndh0004/code/coge_tools/dc_tools/simple_bed_merge.py` to merge adjacent regions with the same 
 call. This is just so we can find The space between regions.
 
 
 ```bash
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ awk 'OFS="\t" {print $2,$1}' ../data/new >cor_genome_sizes.txt
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ head cor_genome_sizes.txt 
Super-Scaffold_1	2385366
Super-Scaffold_10	5195770
Super-Scaffold_100	4772068
Super-Scaffold_1001	1401236
Super-Scaffold_1005	1228904
Super-Scaffold_1007	721261
Super-Scaffold_1008	253901
Super-Scaffold_1009	258217
Super-Scaffold_101	3620639
Super-Scaffold_1015	979394
 ```
  Okay now lets remove ambig from pymerged.
  ```bash
  ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -v ambig pymerged_sorted_xo_chroms_oct1_2018.bed \
   >no_ambig_pymerged_sorted_xo_chroms_oct1_2018.bed
 
 
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep  'ambig' pymerged_sorted_xo_chroms_oct1_2018.bed  \
> ambig_pymerged_sorted_xo_chroms_oct1_2018.bed 
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ bedtools complement  \
-i no_ambig_pymerged_sorted_xo_chroms_oct1_2018.bed -g cor_genome_sizes.txt \
> comp_full_no_ambig_pymerged_sorted_xo_chroms_oct1_2018.bed 

 
 
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f ab_scaf_list.txt  \
 comp_full_no_ambig_pymerged_sorted_xo_chroms_oct1_2018.bed | awk '$2 !=$3'
Super-Scaffold_10	203641	356558
Super-Scaffold_10	2584092	3387643
Super-Scaffold_10	5041777	5195770
Super-Scaffold_1173	518777	1324513
Super-Scaffold_12	3850311	4266075
Super-Scaffold_13	609845	1457381
Super-Scaffold_146	1194073	1950407

 
  ```
  Now we have a set of target sites to validate assemblies in the ambig regions and the spaces between AB and calls.
  We can do a merge and put all these regions together in bedtools so we have a set of target regions to consider. 
  
  ```bash
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ cat ambig_pymerged_sorted_xo_chroms_oct1_2018.bed  \
redu_comp__no_ambig_pymerged_sorted_xo_chroms_oct1_2018.bed |\
 sort -k1 -k2 |\
  awk 'OFS="\t" {print $1,$2,$3}'\
  > total_junction_unmerged_unsorted.bed


ndh0004@IorekByrnison:~/code/coge_tools/data_out$ bedtools sort -i total_junction_unmerged_unsorted.bed > total_junction_unmerged.bed 

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ 

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ wc -l total_junction.bed 
43 total_junction.bed

#43 regions to check for proper assem.

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w Super-Scaffold_76 cor_genome_sizes.txt total_junction.bed 
cor_genome_sizes.txt:Super-Scaffold_76	    6,443,148
total_junction.bed:Super-Scaffold_76	0	6,069,511


```
