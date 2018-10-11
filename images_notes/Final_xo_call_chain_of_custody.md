##Final calls for xo regions.
What we want to do first is to validate the assemblies in the regions of interest.
This can be useful for Hatekayama et al. (2018) data and for Hittalmani et al. (2017).


Analysis was run using 
commits used: b47cc18828ed91f59b041728e90eece40974efeb
            : bd0dc6662b06c74ae0c9de5b537cfc5965c04d1a # for merging mate pair reads later in pipeline.
bedtools v2.25.0
fastalength from exonerate version 2.4.0

bash


First we want to begin by getting a list of all AB contigs. 
I have produced notes on the underlying calls we have made.

[notes_on_AB](https://github.com/NDHall/coge_tools/blob/master/images_notes/process_for_calling_AvsB.md)
[notes_on_dc_ksbagofgenes](https://github.com/NDHall/coge_tools/blob/master/images_notes/notes_dc_ks_bagofgenes.md)
[notes_on_xo_sep24](https://github.com/NDHall/coge_tools/blob/master/images_notes/notes_on_xo_sep24.md)



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
$ #data from new file was previously generenrated on gate using
$ fastalength -v
fastalength from exonerate version 2.4.0
Using glib version 2.48.1
Built on May 18 2016
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

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ git rev-parse HEAD
b47cc18828ed91f59b041728e90eece40974efeb
```
##Mapping
I tried various approaches. The most succesful one was to map uncleaned reads to the entire genome using local and 
take best alignment.  This allows us to look mismatched reads as well. There was some mismatched read support across  

```bash 
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ samtools view -H MP_looseMapLocalNu30K.bam |egrep '^@PG'
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.9	CL:"/automnt/opt/bin/../bowtie2/default/bowtie2-align-s \
--wrapper basic-0 \
--threads 15 \
--local \
--rf \
-X 30000 \
-x Ecor_PR202_scaffolds.fa \
--passthrough \
-1 MP_raw_DRR095902_1.fastq 
-2 MP_raw_DRR095902_2.fastq
```


Next bam was processed with bedtools to create a bed file of mate pair positions and matepair reads were mereged into 
single bed features if they were on the same chromosome so we can view mp support of XO assembly. Score in this version
is based on mapping quality.

```bash 
$   bedtools bamtobed -i MP_looseMapLocalNu30K.bam >MP_looseMapLocalNu30K.bed 
$   /home/ndh0004/code/venv_btools/bin/python3.5 /home/ndh0004/code/coge_tools/viz/bam2bed_mate_merge.py
    Finshed reading Reads.

    5027549 read_matches found

Process finished with exit code 0

$   bedtools sort -i MP_looseMapLocalNu30Kmerged.bed > sorted_MP_looseMapLocalNu30Kmerged.bed 
$ bedtools intersect -a sorted_MP_looseMapLocalNu30Kmerged.bed  -b total_junction.bed -u  >tj_sorted_MP_looseMapLocalNu30Kmerged.bed 
    #-u pulls all bed regions that overlap with juction by as little 1 bp. 
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ git rev-parse HEAD
bd0dc6662b06c74ae0c9de5b537cfc5965c04d1a


```
##Examples:

Green track shows solid connection. which contains multiple high quality mapping matepairs.


###SS95
![Super-Scaffold_95][95_xo]




###SS39 with reads overlapping ambiguous cross over region shown as list. 
![Super-Scaffold_39][39_xo]



[95_xo]:mp_xo_support_95.png
[39_xo]:mp_xo_39.png



##Other visualized regions.
chrom   mp support
* Super-Scaffold_1173 good
* Super-Scaffold_12 good Also interesting overlaps in syntenic regions. These are facilitated by the bag of genes
  approach. Interesting and good that they agree. 
* Super-Scaffold_128 fair. The xo appears to happen inside the first sytenic region which is called as A. Good that
  the mate pairs support A to B call junction but no support at actual cross over. 
  * at the actual crossover it is a little grim.  I did find good support here with tablet viewer. However, There 
  were also several that I moused over in tablet that linked to alt assemblies. This could be gene to gene mapping. But
  I would give this cross over fair. It is one to check perhaps.
* Super-Scaffold_2765 380641 432770  is the linked region to 128 I have changed it to \*\*Provsional after vis.
* Super-Scaffold_13 good
* Super-Scaffold_146 good There is an apparent break in mp coverage
  however this is not inbetween A or B Calls.
* Super-Scaffold_147 good.  
* Super-Scaffold_158 good 
* Super-Scaffold_202 good
* Super-Scaffold_210 crossover inside ambig region
    * region Super-Scaffold_210_2719618_2936984
* Super-Scaffold_2412  cross over inside A called region
    * Super-Scaffold_2412:3849834-4101556
* Super-Scaffold_25 good
* Super-Scaffold_266 good
    * Super-Scaffold_266:1943804-2228668 might contain cross over at
      3' end 
           * ambig Super-Scaffold_266:11783-135556 This is a fine we could even
     call it B by eye. its 4 to 1. Not good enough for chi square
     but plenty strong given its context.
        * also interesting coracana vs coracana blocks are not tracking 
        all the coracana indica all calls are the same over interleaved region
   
* Super-Scaffold_50 
    * cross_over in ambig Super-Scaffold_50:944049-1128041
    * other cross over in inter-syntenic region

* Super-Scaffold_70 mid ambig cross over again
    * Super-Scaffold_70:3321498-3790880
 * Super-Scaffold_76 Cross over in inter-syntentic region
 
* Super-Scaffold_84  This is really interesting with lots of potential 
    * intra-syntentic region cross overs.
      *  Super-Scaffold_84:1082994-124247  
      *  Super-Scaffold_84:2602493-2884299
      *  Super-Scaffold_84:3414116-3478591
      *  Super-Scaffold_84:4040135-4273865
* Super-Scaffold_95 good this is one of the first ones 
    I checked
    * ambig region cross over 
    Super-Scaffold_95:1896657-2247103
* scaffold306_size1284463_obj 
    * interesting cross over. appears pretty sloppy at 
    3' end of AB xo. but clean at 5' Strong mate-pair support.
       it appears there is a problem with the naming convention, pipe vs
       no pipe. 
    * scaffold591_size730508_subseq_1:450970_obj Good 
    * scaffold93_size1230182_subseq_1:987274_obj  crossover mid B call
        * XO scaffold93_size1230182_subseq_1%3A987274_obj:284858..502764 (+ strand)
    




  ![Super-Scaffold_128 XO tablet][xo_128]




[xo_128]:mp_xo_ss128.png 
##Table of Mate Pair Support of XO regions. 



* table does not want to format on preview

|chrom |start|stop|alt|current |trust |
|------|----:|---:|--:|-------:|-----:|
|Super-Scaffold_10|203641|356558|181|480|Trust|
|Super-Scaffold_10|2584092|3387643|807|1837|Trust|
|Super-Scaffold_10 |4944718 |5195770 |390 |741 |Trust|
|Super-Scaffold_105 |3035427 |4238156 |1590 |3751 |Trust|
|Super-Scaffold_1173 |518777 |1324513 |926 |2483 |Trust|
|Super-Scaffold_12 |3850311 |4266075 |547 |1324 |Trust|
|Super-Scaffold_128 |289384 |318912 |21 |88 |\*\*Provisional|
|Super-Scaffold_13 |609845 |1457381 |975 |2337 |Trust|
|Super-Scaffold_146 |1194073 |1950407 |883 |1997 |Trust|
|Super-Scaffold_146 |2263235 |2296778 |31 |130 |Trust|
|Super-Scaffold_147 |1729977 |1904642 |114 |523 |Trust|
|Super-Scaffold_147 |2044366 |2722867 |678 |2100 |Trust|
|Super-Scaffold_158 |393336 |686942 |436 |740 |Trust|
|Super-Scaffold_158 |1117796 |1530254 |436 |1211 |Trust|
|Super-Scaffold_158 |1610047 |1662449 |93 |180 |Trust|
|Super-Scaffold_202 |1644842 |5133839 |2405 |7369 |Trust|
|Super-Scaffold_202 |5186411 |10198086 |3586 |9861 |Trust|
|Super-Scaffold_210 |1999606 |2284541 |214 |579 |Trust|
|Super-Scaffold_210 |2425327 |3210991 |819 |2312 |Trust|
|Super-Scaffold_2412 |3806141 |3849832 |64 |171 |Trust|
|Super-Scaffold_25 |0 |1006350 |466 |1211 |Trust|
|Super-Scaffold_25 |1712239 |1812204 |100 |280 |Trust|
|Super-Scaffold_25 |2131292 |2208811 |54 |235 |Trust|
|Super-Scaffold_253 |2821772 |2983973 |252 |482 |Trust|
|Super-Scaffold_266 |0 |1134602 |1051 |2904 |Trust|
|Super-Scaffold_266 |2228668 |2307550 |38 |250 |Trust|
|Super-Scaffold_2765 |0 |151416 |254 |327 |Trust|
|Super-Scaffold_2765 |380641 |432770 |25 |192 |\*\*Provsional|
|Super-Scaffold_291 |5748370 |8426440 |3862 |8202 |Trust|
|Super-Scaffold_307 |466447 |1181080 |831 |2054 |Trust|
|Super-Scaffold_311 |325086 |1013238 |699 |1718 |Trust|
|Super-Scaffold_39 |2889574 |3609216 |610 |1888 |Trust|
|Super-Scaffold_39 |5389077 |5614697 |228 |625 |Trust|
|Super-Scaffold_39 |5734095 |6249835 |628 |1423 |Trust|
|Super-Scaffold_50 |0 |1286993 |1220 |3177 |Trust|
|Super-Scaffold_50 |5219358 |5804485 |762 |1592 |Trust|
|Super-Scaffold_50 |6795737 |6870923 |126 |238 |Trust|
|Super-Scaffold_70 |3252899 |3790880 |575 |1294 |Trust|
|Super-Scaffold_76 |0 |348153 |423 |855 |Trust|
|Super-Scaffold_76 |556798 |1392247 |940 |2274 |Trust|
|Super-Scaffold_84 |599184 |763955 |75 |948 |Trust|
|Super-Scaffold_84 |950408 |1343812 |260 |1721 |Trust|
|Super-Scaffold_84 |2602492 |2948992 |188 |2739 |Trust|
|Super-Scaffold_84 |3315551 |3509449 |275 |2518 |Trust|
|Super-Scaffold_84 |3967116 |4273865 |122 |2659 |Trust|
|Super-Scaffold_84 |6484426 |6623001 |29 |850 |Trust|
|Super-Scaffold_95 |1564112 |1783658 |260 |571 |Trust|
|Super-Scaffold_95 |1884883 |2384454 |540 |1237 |Trust|
|Super-Scaffold_95 |4880117 |5042869 |152 |368 |Trust|
|scaffold306_size1284463_obj |785065 |964483 |0 |0 |Provisional|
|scaffold591_size730508_subseq_1:450970_obj |208811 |260000 |0 |0 |Provisional|
|scaffold93_size1230182_subseq_1:987274_obj |502766 |536202 |0 |0 |Provisional|
|scaffold93_size1230182_subseq_1:987274_obj |813027 |987274 |0 |0 |Provisional|

 \*\*Downgraded after visual inspections
##The code that I used to look at  cross over support.


```bash 
$ # to get non supporting reads I used the awk statement below.
ndh0004@venus:/home/biolinux/scratch/Eleusine_Synteny_Proj/mapping/Hatakeyama_validation$ \
awk '$1 ~"@" || $7 != "=" ' MP_looseMapLocalNu30K.sam  |\
 samtools view -hb - | \
 samtools sort - MP_looseMapLocalNu30K_alt_chrom &&\
 samtools index MP_looseMapLocalNu30K_alt_chrom.bam & 
$ ~/code/coge_tools/data_out$ bedtools intersect \
-a MP_looseMapLocalNu30K_alt_chrom.bed  \
-b total_junction.bed -u > tj_MP_looseMapLocalNu30K_alt_chrom.bed

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ bedtools coverage -a total_junction.bed  -b tj_sorted_MP_looseMapLocalNu30Kmerged.bed 
Super-Scaffold_10	203641	356558	480	152917	152917	1.0000000
Super-Scaffold_10	2584092	3387643	1837	802492	803551	0.9986821
Super-Scaffold_10	4944718	5195770	741	250809	251052	0.9990321
Super-Scaffold_105	3035427	4238156	3751	1202729	1202729	1.0000000
Super-Scaffold_1173	518777	1324513	2483	805736	805736	1.0000000
Super-Scaffold_12	3850311	4266075	1324	415764	415764	1.0000000
Super-Scaffold_128	289384	318912	88	29528	29528	1.0000000
Super-Scaffold_13	609845	1457381	2337	847536	847536	1.0000000
Super-Scaffold_146	1194073	1950407	1997	756334	756334	1.0000000
Super-Scaffold_146	2263235	2296778	130	33543	33543	1.0000000
Super-Scaffold_147	1729977	1904642	523	174665	174665	1.0000000
Super-Scaffold_147	2044366	2722867	2100	678501	678501	1.0000000
Super-Scaffold_158	393336	686942	740	293606	293606	1.0000000
Super-Scaffold_158	1117796	1530254	1211	412458	412458	1.0000000
Super-Scaffold_158	1610047	1662449	180	52402	52402	1.0000000
Super-Scaffold_202	1644842	5133839	7369	3488997	3488997	1.0000000
Super-Scaffold_202	5186411	10198086	9861	4990531	5011675	0.9957811
Super-Scaffold_210	1999606	2284541	579	284935	284935	1.0000000
Super-Scaffold_210	2425327	3210991	2312	785664	785664	1.0000000
Super-Scaffold_2412	3806141	3849832	171	43691	43691	1.0000000
Super-Scaffold_25	0	1006350	1211	997218	1006350	0.9909256
Super-Scaffold_25	1712239	1812204	280	99965	99965	1.0000000
Super-Scaffold_25	2131292	2208811	235	77519	77519	1.0000000
Super-Scaffold_253	2821772	2983973	482	162201	162201	1.0000000
Super-Scaffold_266	0	1134602	2904	1130548	1134602	0.9964269
Super-Scaffold_266	2228668	2307550	250	78882	78882	1.0000000
Super-Scaffold_2765	0	151416	327	150571	151416	0.9944193
Super-Scaffold_2765	380641	432770	192	52129	52129	1.0000000
Super-Scaffold_291	5748370	8426440	8202	2678070	2678070	1.0000000
Super-Scaffold_307	466447	1181080	2054	714633	714633	1.0000000
Super-Scaffold_311	325086	1013238	1718	688152	688152	1.0000000
Super-Scaffold_39	2889574	3609216	1888	719642	719642	1.0000000
Super-Scaffold_39	5389077	5614697	625	225620	225620	1.0000000
Super-Scaffold_39	5734095	6249835	1423	515740	515740	1.0000000
Super-Scaffold_50	0	1286993	3177	1286441	1286993	0.9995711
Super-Scaffold_50	5219358	5804485	1592	585127	585127	1.0000000
Super-Scaffold_50	6795737	6870923	238	74847	75186	0.9954912
Super-Scaffold_70	3252899	3790880	1294	537981	537981	1.0000000
Super-Scaffold_76	0	348153	855	345535	348153	0.9924803
Super-Scaffold_76	556798	1392247	2274	835449	835449	1.0000000
Super-Scaffold_84	599184	763955	948	164771	164771	1.0000000
Super-Scaffold_84	950408	1343812	1721	393404	393404	1.0000000
Super-Scaffold_84	2602492	2948992	2739	346500	346500	1.0000000
Super-Scaffold_84	3315551	3509449	2518	193898	193898	1.0000000
Super-Scaffold_84	3967116	4273865	2659	306749	306749	1.0000000
Super-Scaffold_84	6484426	6623001	850	138575	138575	1.0000000
Super-Scaffold_95	1564112	1783658	571	219546	219546	1.0000000
Super-Scaffold_95	1884883	2384454	1237	499571	499571	1.0000000
Super-Scaffold_95	4880117	5042869	368	162488	162752	0.9983779
scaffold306_size1284463_obj	785065	964483	0	0	179418	0.0000000
scaffold591_size730508_subseq_1:450970_obj	208811	260000	0	0	51189	0.0000000
scaffold93_size1230182_subseq_1:987274_obj	502766	536202	0	0	33436	0.0000000
scaffold93_size1230182_subseq_1:987274_obj	813027	987274	0	0	174247	0.0000000


ndh0004@IorekByrnison:~/code/coge_tools/data_out$ bedtools coverage -a total_junction.bed  -b tj_MP_looseMapLocalNu30K_alt_chrom.bed 
Super-Scaffold_10	203641	356558	181	7432	152917	0.0486015
Super-Scaffold_10	2584092	3387643	807	34997	803551	0.0435529
Super-Scaffold_10	4944718	5195770	390	16388	251052	0.0652773
Super-Scaffold_105	3035427	4238156	1590	74611	1202729	0.0620348
Super-Scaffold_1173	518777	1324513	926	40680	805736	0.0504880
Super-Scaffold_12	3850311	4266075	547	22481	415764	0.0540715
Super-Scaffold_128	289384	318912	21	945	29528	0.0320035
Super-Scaffold_13	609845	1457381	975	43219	847536	0.0509937
Super-Scaffold_146	1194073	1950407	883	38078	756334	0.0503455
Super-Scaffold_146	2263235	2296778	31	1269	33543	0.0378320
Super-Scaffold_147	1729977	1904642	114	5471	174665	0.0313228
Super-Scaffold_147	2044366	2722867	678	29759	678501	0.0438599
Super-Scaffold_158	393336	686942	436	19013	293606	0.0647568
Super-Scaffold_158	1117796	1530254	436	17386	412458	0.0421522
Super-Scaffold_158	1610047	1662449	93	4256	52402	0.0812183
Super-Scaffold_202	1644842	5133839	2405	133848	3488997	0.0383629
Super-Scaffold_202	5186411	10198086	3586	186206	5011675	0.0371544
Super-Scaffold_210	1999606	2284541	214	10210	284935	0.0358327
Super-Scaffold_210	2425327	3210991	819	38240	785664	0.0486722
Super-Scaffold_2412	3806141	3849832	64	2667	43691	0.0610423
Super-Scaffold_25	0	1006350	466	20055	1006350	0.0199285
Super-Scaffold_25	1712239	1812204	100	4674	99965	0.0467564
Super-Scaffold_25	2131292	2208811	54	2955	77519	0.0381197
Super-Scaffold_253	2821772	2983973	252	11606	162201	0.0715532
Super-Scaffold_266	0	1134602	1051	43798	1134602	0.0386021
Super-Scaffold_266	2228668	2307550	38	1890	78882	0.0239598
Super-Scaffold_2765	0	151416	254	10670	151416	0.0704681
Super-Scaffold_2765	380641	432770	25	1219	52129	0.0233843
Super-Scaffold_291	5748370	8426440	3862	122285	2678070	0.0456616
Super-Scaffold_307	466447	1181080	831	36555	714633	0.0511521
Super-Scaffold_311	325086	1013238	699	29087	688152	0.0422683
Super-Scaffold_39	2889574	3609216	610	26854	719642	0.0373158
Super-Scaffold_39	5389077	5614697	228	10550	225620	0.0467600
Super-Scaffold_39	5734095	6249835	628	24626	515740	0.0477489
Super-Scaffold_50	0	1286993	1220	53443	1286993	0.0415255
Super-Scaffold_50	5219358	5804485	762	32866	585127	0.0561690
Super-Scaffold_50	6795737	6870923	126	5373	75186	0.0714628
Super-Scaffold_70	3252899	3790880	575	25075	537981	0.0466095
Super-Scaffold_76	0	348153	423	18715	348153	0.0537551
Super-Scaffold_76	556798	1392247	940	39743	835449	0.0475708
Super-Scaffold_84	599184	763955	75	2822	164771	0.0171268
Super-Scaffold_84	950408	1343812	260	11782	393404	0.0299489
Super-Scaffold_84	2602492	2948992	188	8260	346500	0.0238384
Super-Scaffold_84	3315551	3509449	275	10400	193898	0.0536364
Super-Scaffold_84	3967116	4273865	122	5638	306749	0.0183798
Super-Scaffold_84	6484426	6623001	29	962	138575	0.0069421
Super-Scaffold_95	1564112	1783658	260	10963	219546	0.0499349
Super-Scaffold_95	1884883	2384454	540	22976	499571	0.0459915
Super-Scaffold_95	4880117	5042869	152	7242	162752	0.0444971
scaffold306_size1284463_obj	785065	964483	0	0	179418	0.0000000
scaffold591_size730508_subseq_1:450970_obj	208811	260000	0	0	51189	0.0000000
scaffold93_size1230182_subseq_1:987274_obj	502766	536202	0	0	33436	0.0000000
scaffold93_size1230182_subseq_1:987274_obj	813027	987274	0	0	174247	0.0000000


Lets look at some of the values athat are zero.
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w 'scaffold591_size730508_subseq_1:450970_obj|scaffold93_size1230182_subseq_1:987274_obj| scaffold306_size1284463_obj'  bog_sep28_qaccall_bag_of_genes.tsv | awk '{print $6} ' |sort -k2  | uniq -c 
     29 scaffold18
      6 scaffold19
     20 scaffold220
     13 scaffold517
     11 scaffold89
   # well these are the indica names let's get the coracana reads tehy are linked to.
   

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f tmp bog_sep28_qaccall_bag_of_genes.tsv | awk '{print $3}' | egrep -v -w 'scaffold591_size730508_subseq_1:450970_obj|scaffold93_size1230182_subseq_1:987274_obj| scaffold306_size1284463_obj'| sort | uniq -c 
      5 Super-Scaffold_1173 # already trusted.
     31 Super-Scaffold_146  # already trusted.
     56 Super-Scaffold_204  # not found
     60 Super-Scaffold_2591 # not found
     24 Super-Scaffold_2978 # not found
     13 Super-Scaffold_376  # not found
     47 Super-Scaffold_70   # already trusted
      7 Super-Scaffold_838  # not found

ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w 'Super-Scaffold_204|\
Super-Scaffold_2591|\
Super-Scaffold_2978|\
Super-Scaffold_376|\
Super-Scaffold_838|\
scaffold591_size730508_subseq_1:450970_obj|\
scaffold93_size1230182_subseq_1:987274_obj|\
scaffold306_size1284463_obj '\
  bog_sep28_qaccalled_region_calls.tsv | sort -k2 
# for ref
#call   chrom                                   start   stop    Agenes   Bgenes   chi square p val        chistat
B	scaffold591_size730508_subseq_1:450970_obj	260001	450970	0	25	5.733031437583875e-07	25.0
A	scaffold591_size730508_subseq_1:450970_obj	81909	208810	11	1	0.003892417122778637	8.333333333333334
B	scaffold93_size1230182_subseq_1:987274_obj	284856	502765	4	15	0.011616891434158384	6.368421052631579
A	scaffold93_size1230182_subseq_1:987274_obj	536203	813027	10	0	0.001565402258002549	10.0
A	Super-Scaffold_204	2205	472083	52	1	2.463009590996734e-12	49.075471698113205
ambig	Super-Scaffold_204	483946	522074	5	1	0.10247043485974942	2.6666666666666665
A	Super-Scaffold_204	709186	802161	12	1	0.002281937253315446	9.307692307692308
ambig	Super-Scaffold_204	962299	1061771	1	0	0.31731050786291115	1.0    # here it appears we have lost a fair \
#                                                                                number of genes to replication.
B	Super-Scaffold_2591	1242970	1333849	1	12	0.002281937253315446	9.307692307692308
B	Super-Scaffold_2591	142640	364699	1	24	4.224909405005696e-06	21.16
B	Super-Scaffold_2591	25251	115197	0	15	0.00010751117672950066	15.0
B	Super-Scaffold_2591	469255	971571	1	52	2.463009590996734e-12	49.075471698113205
ambig	Super-Scaffold_2591	978134	1014496	1	5	0.10247043485974942	2.6666666666666665
A	Super-Scaffold_2978	23769	190020	20	1	3.381272563301577e-05	17.19047619047619
B	Super-Scaffold_376	1094648	1181829	2	10	0.020921335337794035	5.333333333333333
B	Super-Scaffold_376	290396	472982	0	24	9.63357008643095e-07	24.0
B	Super-Scaffold_376	488572	604172	1	12	0.002281937253315446	9.307692307692308
B	Super-Scaffold_376	72999	211256	1	11	0.003892417122778637	8.333333333333334


```
Calls underlying each region  look generally 
solid  in general we see these potential cross over sites, are from unlinked regions. 


 

 
