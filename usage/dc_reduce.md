dc_reduce.py filters syntenic blocks by number of genes. It counts the number of genes reported in each block. This
is important, because not all downstream DAGChainer output files report the number of genes in each block in the 
header. For example the output from quota align does not report gene numbers in the headers of syntenic blocks. 

```bash
usage: dc_reduce.py [-h] -i IN_FILE -o OUT -c CUTOFF

Reduce the syntenic blocks reported by raw number of genes contained in the
syntenic block. Now CoGe synfinder reports genes in block that match. So this
is just a post analysis way to filter syntenic block size.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --infile IN_FILE
                        dc file, tetraploid genome in a column and diploid
                        genome in b column
  -o OUT, --outfile OUT
                        outputfile name. No directory required.
  -c CUTOFF, --cutoff CUTOFF
                        minimum acceptable number of genes.

```