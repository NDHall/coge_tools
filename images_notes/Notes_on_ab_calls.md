
 ###Considering mean_ks values returned per syntenic block.
 
 0.5 appears to be the cutoff value we should use. This is not as terribly at odds with previous emmix predictions. 
 Just here we have fewer values and more sensible bins provided by syntenic block. This is discussed at length in 
 notes_on_emmix_in_R.md (There will be some overlap here.)
 ![mean_ks values histogram for *E. coracana* vs *E. coracana*][mean_hist]
 
 If we take the sum of all the lengths calculated as `stop-start=length`.
```R
> l<- sum(dat$length) # for all calculated values.
> l
[1] 867990121
> l/1000000
[1] 867.9901 # ~868 Mbp for all  CvC syntentic blocks.

> red_dat <- subset(dat,mean_ks<=0.5) # as per histogram
> l<- sum(red_dat$length) ; l/1000000
[1] 559.1658 # ~ 559 Mbp to work with 
``` 
###awk inspection

```bash
 > tail -n +2 filter_sort.tsv | awk '$5 <= 0.5' | wc -l # if ks_mean is less than or equal to new cutoff.
500
> tail -n +2 filter_sort.tsv | awk '$5 <= 0.5 {print $2}' | sort | \
  uniq -c | sort -n  |tail # syn per super scaffold
      8 Super-Scaffold_128
      8 Super-Scaffold_210
      8 Super-Scaffold_266
      8 Super-Scaffold_307
      8 Super-Scaffold_345
      8 Super-Scaffold_514
     10 Super-Scaffold_4946
     11 Super-Scaffold_191
     12 Super-Scaffold_395
     22 Super-Scaffold_152
     
>tail -n +2 filter_sort.tsv | awk '$5 <= 0.5 {print $2}' | sort | uniq -c | sort -n  | wc -l \
# num unique super scafs.
148

``` 
###Bed creation 
Now that we have looked at the data we can go ahead and create a bed of regions to include in our syntentic block 
analysis. Also note that younger syntentic blocks tend to have more genes.



![mean_ks ~ gene_num][ks_by_gene]


There appears to be an outlier. We can pull that one in using `subset()` in R.  

```R
> subset(dat,(num_genes >100) & (mean_ks >= 0.5))
     length              chrom num_genes num_kska_genes mean_ks mean_ka
381 4431440 Super-Scaffold_291       235              1  1.7886  1.6787
382 3422267   Super-Scaffold_8       235              1  1.7886  1.6787
935 3422267   Super-Scaffold_8       235              1  1.7886  1.6787
936 4431440 Super-Scaffold_291       235              1  1.7886  1.6787

#this one exhibits netural selection at a high rate of substitition. It could be a pseudogene or a transposom.
#All 4 entries share the same value to be mapped which mean ks and num_genes. If you look closely you can see
# they are in fact the same regions. This doubling is expected given it is a self vs self comparison.

```

Now we can create a bed file using `awk`

```bash
> awk '$7 <= 0.5 && $7 != "NA" {print $1"\t"$2"\t"$3}' CvC_ks.tsv > ks_include.bed

> wc -l ks_include.bed 
500
```

Let's look at region files produced during this analysis.
```bash

> pwd 
/home/ndh0004/Documents/coracana_AB/full_run_sept10
>ls *region*
ab_calls_CvsI_10cutoff_region_calls.tsv  ab_calls_CvsI_12cutoff_region_calls.tsv  ab_calls_CvsI_14cutoff_region_calls.tsv  ab_calls_CvsI_5cutoff_region_calls.tsv  ab_calls_CvsI_8cutoff_region_calls.tsv

> cat *region* |awk '{print $2,$3,$4}'| sort | uniq -c | head 
      2 scaffold1083_size232329_obj 26 220961
      1 scaffold112_size2445781_subseq_1:195369_obj 50904 147964
      1 scaffold1208_size198491_obj 95584 165033
      1 scaffold1379_size156809_obj 7500 62050
      5 scaffold1513_size128307_obj 22 97514
      2 scaffold1517_size127354_obj 4741 109542
      3 scaffold1536_size123903_obj 459 97157
      3 scaffold1599_size113892_obj 22183 103488
      1 scaffold1605_size112729_obj 36216 82095
      2 scaffold1719_size93680_obj 21235 92994


>  cat *region* |awk '{print $2,$3,$4}'| sort | uniq -c | awk '{print $1}' | sort | uniq -c  
#print out and count unique bed regions.
    510 1
    234 2
    120 3
     88 4
     64 5


> cat *region* |awk '{print $2}'| sort | uniq -c | sort -n |tail 
     37 Super-Scaffold_10
     37 Super-Scaffold_311
     38 Super-Scaffold_15
     38 Super-Scaffold_19
     38 Super-Scaffold_210
     38 Super-Scaffold_84
     45 Super-Scaffold_24
     46 Super-Scaffold_12
     55 Super-Scaffold_50
     65 Super-Scaffold_18

# in the next 2 examples you can verify by eye that we are not pulling in any overlaps. This is 
# a function of the the script.

> cat *region* |awk '{print $2,$3,$4}'| sort | uniq -c | egrep -w 'Super-Scaffold_8'
      1 Super-Scaffold_8 1373958 1422464
      1 Super-Scaffold_8 2414272 2514231
      2 Super-Scaffold_8 2680151 2782102
      1 Super-Scaffold_8 2920534 3035181
      1 Super-Scaffold_8 394682 435285
      1 Super-Scaffold_8 500348 660293
      2 Super-Scaffold_8 777791 1221206
      
 >egrep 'Super-Scaffold_8.*1373958.*1422464' *region* # check which file a uniq hit occurs in.
ab_calls_CvsI_5cutoff_region_calls.tsv:B	Super-Scaffold_8	\
1373958	1422464	0	5	0.025347318677468325	5.0


> cat *region* |awk '{print $2,$3,$4}'| sort | uniq -c | egrep -w 'Super-Scaffold_18'
      2 Super-Scaffold_18 1194533 1238316
      2 Super-Scaffold_18 1252344 1368670
      1 Super-Scaffold_18 1420649 1466922
      1 Super-Scaffold_18 1621764 1820015
      5 Super-Scaffold_18 1842562 1993197
      1 Super-Scaffold_18 2073796 2147232
      5 Super-Scaffold_18 2189498 2331558
      1 Super-Scaffold_18 2539471 2622753
      1 Super-Scaffold_18 2675925 2774743
      4 Super-Scaffold_18 2798672 3015957
      4 Super-Scaffold_18 310431 454443
      3 Super-Scaffold_18 3126222 3375044
      3 Super-Scaffold_18 3385174 3457324
      2 Super-Scaffold_18 3499937 3568756
      2 Super-Scaffold_18 3580534 3640681
      1 Super-Scaffold_18 3661879 3851503
      1 Super-Scaffold_18 4026870 4145043
      1 Super-Scaffold_18 4325936 4413439
      2 Super-Scaffold_18 4415669 4621308
      4 Super-Scaffold_18 4704292 4911819
      1 Super-Scaffold_18 4912938 5052859
      1 Super-Scaffold_18 6088808 6191043
      1 Super-Scaffold_18 6434384 6451953
      3 Super-Scaffold_18 6478143 6641908
      1 Super-Scaffold_18 6678064 6706232
      4 Super-Scaffold_18 6755043 6999747
      3 Super-Scaffold_18 7050831 7128643
      2 Super-Scaffold_18 7140866 7167691
      3 Super-Scaffold_18 87018 306416

> cat *region* |awk 'BEGIN{OFS="\t";} {print $2,$3,$4}'| sort | uniq  | wc -l 
1016

> cat *region* |awk 'BEGIN{OFS="\t";} {print $2,$3,$4}'| sort | uniq  > total_regions.bed

> qsub -I -V 
[ndh0004@node124 full_run_sept10]$ source ~/Documents/coracana_AB/venv/bin/activate
# now in venv
(venv) [ndh0004@node124 full_run_sept10]$ python ~/src/coge_tools/dc_tools/strict_regions_to_chrom_level_calls.py  -h 
usage: strict_regions_to_chrom_level_calls.py [-h] -i IN_FILE -o OUT [-k KEEP]
                                              [-e EXCLUDE]

This script takes a series of region calls and filters to include targeted bed
regions or exclude targeted bed regions.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --infile IN_FILE
                        region_calls file produced by ab_call.py
  -o OUT, --outfile OUT
  -k KEEP, --regions_to_keep KEEP
  -e EXCLUDE, --regions_to_exclude EXCLUDE


# needs region file not a bed so I uniqued all region files. 
# I got this result that needs to be looked at. We have already established that
# regions are correctly called above. But the whole file was giving me a weird case
#where it was off by one line on the uniq between reduced bed and the whole file. 
```

####A deeper look at problem region.


```bash
 (venv) [ndh0004@node124 full_run_sept10]$ awk '{print $2,$3,$4}' uniq_regions.tsv| sort | uniq -c | sort -n |  tail 
      1 Super-Scaffold_9 744045 832956
      1 Super-Scaffold_980 153303 329951
      1 Super-Scaffold_98 1097382 1443503
      1 Super-Scaffold_99 120402 254719
      1 Super-Scaffold_99 1979675 2034410
      1 Super-Scaffold_99 3053488 3308129
      1 Super-Scaffold_99 3441104 3498518
      1 Super-Scaffold_99 3551256 3600396
      1 Super-Scaffold_9 979861 1090522
      2 Super-Scaffold_179 1461128 1869363
>less -S + search 
A       Super-Scaffold_179      1461128 1869363 29      2       1.238709698569697e-06   23.516129032258064
A       Super-Scaffold_179      1461128 1869363 9       2       0.03480847881186725     4.454545454545454

#gene numbers are off. This needs to be looked at. It appears to be directly related to the number of genes that 
# available for 1 to 1 comparison. It increases as stringency increases and more genes are freed up.
# what is it doing in the underlying network?

 grep scaffold38_ *net| \
 sort -k2 |\
 sed 's/51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.//' \
 | sed 's/[: ]/\t/g;' 
gt10.parse.net	scaffold38_29751_162685	scaffold1517_size127354_obj_4741_109542
gt12.parse.net	scaffold38_29751_162685	scaffold1517_size127354_obj_4741_109542
gt5.parse.net	scaffold38_29751_162685	scaffold1517_size127354_obj_4741_109542
gt8.parse.net	scaffold38_29751_162685	scaffold1517_size127354_obj_4741_109542
gt10.parse.net	scaffold38_220216_491059	Super-Scaffold_1090_41_227049
gt12.parse.net	scaffold38_220216_491059	Super-Scaffold_1090_41_227049
gt14.parse.net	scaffold38_220216_491059	Super-Scaffold_1090_41_227049
gt5.parse.net	scaffold38_220216_491059	Super-Scaffold_1090_41_227049
gt8.parse.net	scaffold38_220216_491059	Super-Scaffold_1090_41_227049
gt10.parse.net	scaffold38_29751_523588	Super-Scaffold_179_1461128_1869363
gt12.parse.net	scaffold38_29751_523588	Super-Scaffold_179_1461128_1869363
gt14.parse.net	scaffold38_29751_523588	Super-Scaffold_179_1461128_1869363
gt5.parse.net	scaffold38_29751_523588	Super-Scaffold_179_1461128_1869363
gt8.parse.net	scaffold38_29751_523588	Super-Scaffold_179_1461128_1869363
gt5.parse.net	scaffold38_153587_281083	Super-Scaffold_24_11411448_11689082
gt5.parse.net	scaffold38_156743_436786	Super-Scaffold_2761_25074_485336
gt8.parse.net	scaffold38_156743_436786	Super-Scaffold_2761_25074_485336


```
###What are the genes doing since they seem to be variable

```bash 
>  grep Super-Scaffold_179_1461128_1869363 \
gene_asyn_aper_bsyn_bper_51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.gt*|\
sort -k2 |\
sed 's/^.*51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.//'\
| sed 's/[: ]/\t/g;'| sort -k2 |\
sed 's/^.*51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.//' |\
 sed 's/[: ]/\t/g;'| cut -b 3-| sed 's/.parse/ parse/' | sort -n
 
 10 parse.txt	g4747	Super-Scaffold_179_1461128_1869363	[99.7]	scaffold1517_size127354_obj_4741_109542	[98.66]
10 parse.txt	g4751	Super-Scaffold_179_1461128_1869363	[99.42]	scaffold1517_size127354_obj_4741_109542	[86.5]
10 parse.txt	g4753	scaffold1517_size127354_obj_4741_109542	[96.03]	Super-Scaffold_179_1461128_1869363	[95.19]
10 parse.txt	g4754	Super-Scaffold_179_1461128_1869363	[98.64]	scaffold1517_size127354_obj_4741_109542	[98.49]
10 parse.txt	g4755	Super-Scaffold_179_1461128_1869363	[100.0]	scaffold1517_size127354_obj_4741_109542	[88.62]
10 parse.txt	g4760	scaffold1517_size127354_obj_4741_109542	[96.55]	Super-Scaffold_179_1461128_1869363	[95.48]
10 parse.txt	g4765	Super-Scaffold_179_1461128_1869363	[97.33]	scaffold1517_size127354_obj_4741_109542	[95.66]
10 parse.txt	g4766	Super-Scaffold_179_1461128_1869363	[99.94]	scaffold1517_size127354_obj_4741_109542	[99.3]
10 parse.txt	g4767	Super-Scaffold_179_1461128_1869363	[100.0]	scaffold1517_size127354_obj_4741_109542	[98.84]
10 parse.txt	g4768	Super-Scaffold_179_1461128_1869363	[99.58]	scaffold1517_size127354_obj_4741_109542	[89.71]
10 parse.txt	g4769	Super-Scaffold_179_1461128_1869363	[100.0]	scaffold1517_size127354_obj_4741_109542	[98.17]
12 parse.txt	g4747	Super-Scaffold_179_1461128_1869363	[99.7]	scaffold1517_size127354_obj_4741_109542	[98.66]
12 parse.txt	g4751	Super-Scaffold_179_1461128_1869363	[99.42]	scaffold1517_size127354_obj_4741_109542	[86.5]
12 parse.txt	g4753	scaffold1517_size127354_obj_4741_109542	[96.03]	Super-Scaffold_179_1461128_1869363	[95.19]
12 parse.txt	g4754	Super-Scaffold_179_1461128_1869363	[98.64]	scaffold1517_size127354_obj_4741_109542	[98.49]
12 parse.txt	g4755	Super-Scaffold_179_1461128_1869363	[100.0]	scaffold1517_size127354_obj_4741_109542	[88.62]
12 parse.txt	g4760	scaffold1517_size127354_obj_4741_109542	[96.55]	Super-Scaffold_179_1461128_1869363	[95.48]
12 parse.txt	g4765	Super-Scaffold_179_1461128_1869363	[97.33]	scaffold1517_size127354_obj_4741_109542	[95.66]
12 parse.txt	g4766	Super-Scaffold_179_1461128_1869363	[99.94]	scaffold1517_size127354_obj_4741_109542	[99.3]
12 parse.txt	g4767	Super-Scaffold_179_1461128_1869363	[100.0]	scaffold1517_size127354_obj_4741_109542	[98.84]
12 parse.txt	g4768	Super-Scaffold_179_1461128_1869363	[99.58]	scaffold1517_size127354_obj_4741_109542	[89.71]
12 parse.txt	g4769	Super-Scaffold_179_1461128_1869363	[100.0]	scaffold1517_size127354_obj_4741_109542	[98.17]
14 parse.txt	g4783	Super-Scaffold_179_1461128_1869363	[98.83]	Super-Scaffold_1090_41_227049	[97.2]
14 parse.txt	g4784	Super-Scaffold_179_1461128_1869363	[99.26]	Super-Scaffold_1090_41_227049	[97.91]
14 parse.txt	g4785	Super-Scaffold_179_1461128_1869363	[99.8]	Super-Scaffold_1090_41_227049	[98.26]
14 parse.txt	g4787	Super-Scaffold_179_1461128_1869363	[95.13]	Super-Scaffold_1090_41_227049	[88.98]
14 parse.txt	g4788	Super-Scaffold_179_1461128_1869363	[99.91]	Super-Scaffold_1090_41_227049	[88.2]
14 parse.txt	g4790	Super-Scaffold_179_1461128_1869363	[96.15]	Super-Scaffold_1090_41_227049	[86.75]
14 parse.txt	g4791	Super-Scaffold_179_1461128_1869363	[99.75]	Super-Scaffold_1090_41_227049	[83.59]
14 parse.txt	g4792	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[95.83]
14 parse.txt	g4793	Super-Scaffold_179_1461128_1869363	[99.75]	Super-Scaffold_1090_41_227049	[97.84]
14 parse.txt	g4794	Super-Scaffold_179_1461128_1869363	[99.49]	Super-Scaffold_1090_41_227049	[98.79]
14 parse.txt	g4795	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[97.96]
14 parse.txt	g4796	Super-Scaffold_179_1461128_1869363	[99.66]	Super-Scaffold_1090_41_227049	[96.28]
14 parse.txt	g4797	Super-Scaffold_179_1461128_1869363	[96.53]	Super-Scaffold_1090_41_227049	[94.61]
14 parse.txt	g4798	Super-Scaffold_179_1461128_1869363	[99.43]	Super-Scaffold_1090_41_227049	[98.47]
14 parse.txt	g4799	Super-Scaffold_179_1461128_1869363	[99.2]	Super-Scaffold_1090_41_227049	[99.11]
14 parse.txt	g4801	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[98.49]
14 parse.txt	g4803	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[91.43]
14 parse.txt	g4807	Super-Scaffold_179_1461128_1869363	[98.06]	Super-Scaffold_1090_41_227049	[94.89]
14 parse.txt	g4808	Super-Scaffold_179_1461128_1869363	[99.56]	Super-Scaffold_1090_41_227049	[94.73]
14 parse.txt	g4809	Super-Scaffold_179_1461128_1869363	[99.36]	Super-Scaffold_1090_41_227049	[99.25]
14 parse.txt	g4810	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[97.17]
14 parse.txt	g4816	Super-Scaffold_179_1461128_1869363	[98.02]	Super-Scaffold_1090_41_227049	[93.64]
14 parse.txt	g4817	Super-Scaffold_179_1461128_1869363	[99.74]	Super-Scaffold_1090_41_227049	[97.8]
14 parse.txt	g4818	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[98.1]
14 parse.txt	g4819	Super-Scaffold_179_1461128_1869363	[99.51]	Super-Scaffold_1090_41_227049	[96.01]
14 parse.txt	g4820	Super-Scaffold_1090_41_227049	[96.34]	Super-Scaffold_179_1461128_1869363	[95.01]
14 parse.txt	g4821	Super-Scaffold_179_1461128_1869363	[87.06]	Super-Scaffold_1090_41_227049	[85.12]
14 parse.txt	g4822	Super-Scaffold_179_1461128_1869363	[99.28]	Super-Scaffold_1090_41_227049	[98.26]
14 parse.txt	g4823	Super-Scaffold_1090_41_227049	[98.63]	Super-Scaffold_179_1461128_1869363	[98.05]
14 parse.txt	g4829	Super-Scaffold_179_1461128_1869363	[99.54]	Super-Scaffold_1090_41_227049	[93.91]
14 parse.txt	g4830	Super-Scaffold_179_1461128_1869363	[100.0]	Super-Scaffold_1090_41_227049	[95.12]


#Let's look at what genes are changing among calls. 
>grep Super-Scaffold_179_1461128_1869363 \
gene_asyn_aper_bsyn_bper_51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.gt*|\
sort -k2 |\
sed 's/^.*51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.//'|\
 sed 's/[: ]/\t/g;'|\
 sort -k2 |\
 sed 's/[: ]/\t/g;'|\
 cut -b 3-|\
 sed 's/.parse/ parse/' |\ sort -n   \
| awk '{print $3}' |\
 sort | uniq -c | sort -n 
      1 g4783
      1 g4784
      1 g4785
      1 g4787
      1 g4788
      1 g4790
      1 g4791
      1 g4792
      1 g4793
      1 g4794
      1 g4795
      1 g4796
      1 g4797
      1 g4798
      1 g4799
      1 g4801
      1 g4803
      1 g4807
      1 g4808
      1 g4809
      1 g4810
      1 g4816
      1 g4817
      1 g4818
      1 g4819
      1 g4820
      1 g4821
      1 g4822
      1 g4823
      1 g4829
      1 g4830
      2 g4747
      2 g4751
      2 g4753
      2 g4754
      2 g4755
      2 g4760
      2 g4765
      2 g4766
      2 g4767
      2 g4768
      2 g4769


```
Now create the region file for calls to be made from. 


[raw_hist]:/home/ndh0004/code/coge_tools/images_notes/raw_hist.png
[lg_hist]:/home/ndh0004/code/coge_tools/images_notes/focused_hist_lg.png
[dn_plt]:/home/ndh0004/code/coge_tools/images_notes/dens_vline.png
[loc_min]:/home/ndh0004/code/coge_tools/images_notes/local_min.png
[mean_hist]:mean_ks_vals.png
[ivc]:1.png
[cvc]:cvc_coge.png
[ks_by_gene]:mks_ng.png