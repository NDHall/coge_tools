###dc_gene_count.py 

This script produces a fairly basic histogram of genes per syntentic block. It can be helpful, because even though 
you establish a minimum cutoff in CoGe it is not always apparent how gene numbers from syntenic blocks are distrbuted. 
We found this helpful in determining basic cutoffs to eliminate unwanted WGD signatures.

```bash
Takes a DAG_Chainer(dc) file and returns a historgram of number of genes in
each syntenic block.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --infile IN_FILE
                        dc file
  -d OUT_DIR, --dir OUT_DIR
                        directory to save results to, default is current
                        working directory
  --dpi DPI             dpi for output graph (default == 300)
  --bins BINS           size of bins for histogram (default == 20)
  -s, --show            show graph when made (default is False)
  -o OUT, --outfile OUT

> python dc_tools/dc_gene_count.py -i test_data/SScaf_indica_v_cor.test \
                                   -o test_out/sscaf_testv5

```

![ output from gene _count.py run with default settings ][out]


[out]:sscaf_testv5.png "output from gene _count.py run with default settings"


