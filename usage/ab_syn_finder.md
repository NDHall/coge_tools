*ab_syn_finder.py* takes DAGChainer output from CoGe SynMap and finds non-overlapping syntenic blocks that occur in 
a 2:1 (allopolypoid:parent) ratio. The filtering used in this script is stricter than that used by Quota Align since
it does not permit syntenic blocks to overlap on the sequence. This is addressed in more detail in 
*notes_on_qac_aug30.md*. If you would like to rely heavily on Quota Align fitering, you should also look at the 
stringency that is being applied to calling of syntenic blocks. It should be possible to adjust the stringency such 
that overlapping regions are minimized. This is also the best practice for *ab_syn_finder.py* since it throws away all
overlapping regions. It is possible to lose A vs B calls because the intial syntenic block calls were too permisive. 

```bash
usage: ab_syn_finder.py [-h] -i IN_FILE [-d OUT_DIR] -o OUT_PRE

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --infile IN_FILE
                        dc file, tetraploid genome in a column and diploid
                        genome in b column
  -d OUT_DIR, --dir OUT_DIR
                        directory to save results to, default is current
                        working directory
  -o OUT_PRE, --outfile OUT_PRE

```
