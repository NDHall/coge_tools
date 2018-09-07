if we cat all syn_ab_call files together. We get interesting results. Each should be have adhered to the same guidelines
but the behavior of syntenic block retention will differ because syntenic block purity will be affected differentially
at each new cutoff. 

```bash
[ndh0004@hopper-login](full_run_sept6)[16:26]: wc -l *_all*tsv |sort -n 
  136 ab_calls_CvsI_12cutoff_chrom_calls_all.tsv
  141 ab_calls_CvsI_10cutoff_chrom_calls_all.tsv
  167 ab_calls_CvsI_8cutoff_chrom_calls_all.tsv
  202 ab_calls_CvsI_5cutoff_chrom_calls_all.tsv
  646 total

```
Even though more calls were made with 5 cutoff how many were good?

```bash
[ndh0004@hopper-login](full_run_sept6)[16:26]: wc -l *_good*tsv |sort -n 
  106 ab_calls_CvsI_10cutoff_chrom_calls_good.tsv
  106 ab_calls_CvsI_12cutoff_chrom_calls_good.tsv
  119 ab_calls_CvsI_8cutoff_chrom_calls_good.tsv
  130 ab_calls_CvsI_5cutoff_chrom_calls_good.tsv
  461 total

```

are they all the same? 

```bash
[ndh0004@hopper-login](full_run_sept6)[16:35]: cat *_good*tsv |\
 awk '{print $5}' | sort | uniq | sort -n  | wc -l 
217

[ndh0004@hopper-login](full_run_sept6)[16:35]: cat *_good*tsv |\
 awk '{print $1,$5}' | sort | uniq | sort -n  | wc -l 
224

```
No they are not always the same 7 have different calls among all the different cutoffs. However, the vast majority,
when they appear are called identicaly.


 ```bash

[ndh0004@hopper-login](full_run_sept6)[16:32]: cat *_good*tsv |\
 awk '{print $1,$5}' | sort | uniq | awk '{print $2}' | sort | uniq -c | sort -n |tail  
      1 Super-Scaffold_98
      1 Super-Scaffold_980
      1 Super-Scaffold_99
      2 Super-Scaffold_10
      2 Super-Scaffold_100
      2 Super-Scaffold_266
      2 Super-Scaffold_291
      2 Super-Scaffold_307
      2 Super-Scaffold_311
      2 Super-Scaffold_95

```
Most calls remain the same however we see that a small subset of them actually are not the same across all sets.
This does suggest that results are fairly robust, but not always perfect. It is possible then to take all these calls 
together if we are not double dipping in terms of gene analysis. That is every A vs B call is made once and only per 
gene. And these genes do not float between regions, based on stringency. There are a couple of nice features to the 
coding that allow us to take all these together. First, overlap calls are made on a completed DAGChainer bed object. 
This means that every region considered in all files is not overlapping, and occurs at a 2:1 ratio. The reason we 
recover more unique regions as we increase stringency is that we slowly eliminate shorter syntenic block calls which 
are, in plants likely relictual signatures of ancient, and not so ancient WGD events. This is expected. It also means 
the shorter nested in the presumably larger and more recent blocks are never counted. As we raise the minimum, specious
hits are weeded out and we recover larger blocks. 
 
![alt text]()



