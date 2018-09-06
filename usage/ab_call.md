
###inputs
ab_call.py takes as output a series of non-overlapping syntenic blocks. These blocks must occur at a t 2:1 ration of 
allopolyploid to parent. Sample input file where col1 == region id , col2 == AvsB, col3 == geneid. 

```bash
Super-Scaffold_50_1286993_1459256       A       g6833
Super-Scaffold_12_100830_248699 B       g6833
Super-Scaffold_50_1286993_1459256       A       g6827
Super-Scaffold_12_100830_248699 B       g6827
Super-Scaffold_50_1286993_1459256       A       g6828
Super-Scaffold_12_100830_248699 B       g6828
Super-Scaffold_50_1286993_1459256       A       g6844
Super-Scaffold_12_100830_248699 B       g6844
Super-Scaffold_50_1286993_1459256       A       g6837
Super-Scaffold_12_100830_248699 B       g6837

```
###calls
ab_call.py parses this file and tallies A vs B by region and then by chrom. Final chrom calls are informed by the 
strength of A vs B calls across entire chrom. 
```bash
usage: ab_call.py [-h] -i IN_FILE -o OUT_PRE [-p PMIN] [-c CUTOFF]

takes syn_block output from ab_synfinder.py and calls A vs B for each region.
If regions fall outside of cutoffs for regions output ambig is designated. If
this occurs at the chrom/contig level, ambiguous calls can be found in
_ambig.tsv

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --infile IN_FILE
                        per gene call from ab_syn_finder.py
  -o OUT_PRE, --outfile OUT_PRE
                        ab_call.py produces mulitple out files each with their
                        own name. Here you provided the base of the name
                        ab_call.py will use to write out files for this
                        analysis.
  -p PMIN, --pval_cutoff PMIN
                        Minimum acceptable p-value for a per region call. In
                        which chi square analysis is used to determine A vs B
                        call. default == 0.05
  -c CUTOFF, --cutoff CUTOFF
                        Minimum acceptable number of genes used to call AvsB.
                        It is suggested that this value is informed by
                        examination of WGD signatures from parent genomes.
                        default == 15

```

###outputs
For example an -o test_v4a will produce the following files. Notice reg_calls.tsv does not have labeled columns. 

```bash
syn_block_ab_gene_test_v4a_chrom_calls_all.tsv
#call   call_qual       num_good_reg    num_ambig_reg   chrom   A_gene_calls    B_gene_calls    frac_A  frac_B  ambig_a ambig_b
AB      strong  10      0       Super-Scaffold_12       37      88      0.296   0.704   0       0
B       strong  2       0       Super-Scaffold_149      2       40      0.047619047619047616    0.9523809523809523      0       0
AB      strong  10      0       Super-Scaffold_50       88      37      0.704   0.296   0       0
A       strong  2       0       Super-Scaffold_37       40      2       0.9523809523809523      0.047619047619047616    0       0
B       strong  2       0       Super-Scaffold_133      1       19      0.05    0.95    0       0
A       strong  2       0       Super-Scaffold_10       19      1       0.95    0.05    0       0
B       provisional     1       1       Super-Scaffold_8        0       5       0.0     1.0     1       4
A       provisional     1       1       Super-Scaffold_291      5       0       1.0     0.0     4       1


syn_block_ab_gene_test_v4a_chrom_calls_ambig.tsv
#call   call_qual       num_good_reg    num_ambig_reg   chrom   A_gene_calls    B_gene_calls    frac_A  frac_B  ambig_a ambig_b
syn_block_ab_gene_test_v4a_chrom_calls_conf.tsv
B       Super-Scaffold_149      2       40      0.047619047619047616    0.9523809523809523
B       Super-Scaffold_133      1       19      0.05    0.95
A       Super-Scaffold_291      9       1       0.9     0.1
B       Super-Scaffold_8        1       9       0.1     0.9
A       Super-Scaffold_10       19      1       0.95    0.05
A       Super-Scaffold_37       40      2       0.9523809523809523      0.047619047619047616


syn_block_ab_gene_test_v4a_chrom_calls_good.tsv
#call   call_qual       num_good_reg    num_ambig_reg   chrom   A_gene_calls    B_gene_calls    frac_A  frac_B  ambig_a ambig_b
AB      strong  10      0       Super-Scaffold_12       37      88      0.296   0.704   0       0
B       strong  2       0       Super-Scaffold_149      2       40      0.047619047619047616    0.9523809523809523      0       0
AB      strong  10      0       Super-Scaffold_50       88      37      0.704   0.296   0       0
A       strong  2       0       Super-Scaffold_37       40      2       0.9523809523809523      0.047619047619047616    0       0
B       strong  2       0       Super-Scaffold_133      1       19      0.05    0.95    0       0
A       strong  2       0       Super-Scaffold_10       19      1       0.95    0.05    0       0


syn_block_ab_gene_test_v4a_chrom_calls_prov.tsv
#call   call_qual       num_good_reg    num_ambig_reg   chrom   A_gene_calls    B_gene_calls    frac_A  frac_B  ambig_a ambig_b
B       provisional     1       1       Super-Scaffold_8        0       5       0.0     1.0     1       4
A       provisional     1       1       Super-Scaffold_291      5       0       1.0     0.0     4       1


syn_block_ab_gene_test_v4a_region_calls.tsv
B       Super-Scaffold_149      803775  1016989 1       21      2.007865612426481e-05   18.181818181818183
B       Super-Scaffold_50       6585961 6795737 2       19      0.0002075015940618336   13.761904761904763
B       Super-Scaffold_8        1373958 1422464 0       5       0.025347318677468325    5.0
A       Super-Scaffold_37       2737691 2966196 21      1       2.007865612426481e-05   18.181818181818183
A       Super-Scaffold_50       3026327 3270181 27      0       2.0345546145444302e-07  27.0
B       Super-Scaffold_12       100830  248699  0       14      0.00018281063298183515  14.0
B       Super-Scaffold_149      440622  654736  1       19      5.699411623331848e-05   16.2
A       Super-Scaffold_12       4977287 5072854 11      0       0.0009111188771537126   11.0
B       Super-Scaffold_50       6808538 6862082 1       7       0.033894853524689295    4.5
B       Super-Scaffold_12       1014303 1170151 0       15      0.00010751117672950066  15.0


```
regions output call chrom start stop A-genes B-genes p-value chi square statistic 


