###Bag of Genes Approach
1. use stringent qac filtering.
2. take ks file and parse into gene calls.
3. It is possible to use either ks with percent id or just percent id.
4. When we call `strict_ks = True` vs `strict_ks = False`

```bash
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ wc -l bog_sep18_qac*_abcalls.tsv 
 12796 bog_sep18_qac_abcalls.tsv
  3078 bog_sep18_qac_strict_abcalls.tsv
 15874 total

```

####Investigation of call files.
```bash
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ awk '{print $2}' bog_sep18_qac_bag_of_genes.tsv | sort | uniq -c | awk '{print $1}' | sort | uniq -c | sort -n 
     25 4
     83 3
   4042 1
   6503 2

```
This suggests we should have 7,006 called anchors
based on number alone. But if we filter by percent id...
```bash
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ awk '$10 >=90 {print $2}' bog_sep18_qac_bag_of_genes.tsv | sort |uniq -c |awk '{print $1}' |sort| uniq -c | sort -n 
     13 4
     42 3
   4107 1
   6398 2
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ wc -l bog_sep18_qac_abcalls.tsv 
12796 bog_sep18_qac_abcalls.tsv
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ echo $((6398 * 2 ))
12796
# again with ks_strict values we can reproduce the number of lines in ab_calls
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ awk '$12 <=3.0 && $10 >=90  {print $2}' bog_sep18_qac_bag_of_genes.tsv | sort |uniq -c |awk '{print $1}' | sort | uniq -c | sort -n 
      6 4
     25 3
   1815 2
   4734 1
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ echo $(( 1815 *2 ))
3630
#3 was not the ks_cutoff
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ awk '$12 <=2.0 && $10 >=90  {print $2}' bog_sep18_qac_bag_of_genes.tsv | sort |uniq -c |awk '{print $1}' | sort | uniq -c | sort -n 
      5 4
     22 3
   1539 2
   4621 1
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ echo $(( 1539 *2 ))
3078
#2 was set as cutoff...

```

```python
    pid_cutoff = 90.0
    ks_cutoff = 2
    syn_len_cutoff = 5
    strict_ks = False
    qac = True
```
The script is behaving as expected for making ab_calls. Notice that we did not filter by
syntenic length `syn_len_cutoff = 5` which is the minimum value reported by CoGe SynMap
If we had chose different lengths this would affect ouf final number of ab_calls.
As it stands now we are letting gene number and qac do the sorting for us. 


You may have also noticed that `qac = True`. a qac.ks file has different section delimiters 
than a normal gcoords.ks file. Setting `qac` to `True` deals with this problem.
 
We also have not lost any lines between qac.ks output and bag_of_genes.tsv 
```bash
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ egrep -v '^#' ../data/51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.Dm0.ma1.qac2.1.50.gcoords.ks | wc -l 
17397
```
and this screen shot from Jupyter
![screen_shot](jupyter_qac_1.png) 

I have now updated 
`parse_ks(infile,pid_cutoff,ks_cutoff,syn_len_cutoff,out_file,call=call,strict_ks=strict_ks, qac=qac)`
to include variable call. Call provides a default label in the 1st column
this will allow us to label multiple files for concatenation. At some 
point adding a seperate file column might be useful. But at this point, we will 
avoid doing that. 
Call update was successfull as suggested by the fact the number of calls don't
change just their labels.  The diffs of whole files change because the 
I changed parameters for calculating `meanks`. 

###Overlap investigation

We want to look at what predicted recombinant regions look like. Some are fairly straight forwards such as 
![ideal][ideal]
There are more complex overlaps such as the next one we will investigate here. 
Let's make sure we have made correct calls. 
![complex][complex]
let's begin by considering what is happening with p_scaffold391

![p_scaffold391][complex_sub1]

The connection between p_scaffold391 and A_SS128 and B_SS128. is a connection made by seven genes. these genes were not
part of the call for A.  but they likely show that exact point of recombination resides in p_scaffold391. 


```bash
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w Super-Scaffold_128 51576_52024_qac_ks_sep19a_calls_5_region_calls.tsv | sort A	Super-Scaffold_128	36603	289383	17	7	0.041226833337163815	4.166666666666667
B	Super-Scaffold_128	2232219	2614392	3	42	6.107886359709498e-09	33.8
B	Super-Scaffold_128	2609494	3115952	1	36	8.717444154664979e-09	33.108108108108105
B	Super-Scaffold_128	3106960	3508023	0	29	7.23782987174e-08	29.0
B	Super-Scaffold_128	318913	524163	1	23	7.09790833090828e-06	20.166666666666668
B	Super-Scaffold_128	3518031	3606981	0	14	0.00018281063298183515	14.0
B	Super-Scaffold_128	3608587	4019641	0	42	9.127341799289188e-11	42.0
B	Super-Scaffold_128	4132382	4273583	1	18	9.616588268883469e-05	15.210526315789474
B	Super-Scaffold_128	4430066	4593260	3	15	0.004677734981047276	8.0
B	Super-Scaffold_128	4883608	4970990	0	12	0.0005320055051392492	12.0
B	Super-Scaffold_128	542218	757393	0	22	2.726504656155499e-06	22.0
B	Super-Scaffold_128	884540	2178431	1	37	5.220985924053899e-09	34.10526315789474
ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w 'scaffold391' 51576_52024_qac_ks_sep20_bag_of_genes.tsv| awk '$1=="A" || $1=="B" {print $1, $14,$15,$16,$17,$18,$19,$3,$4,$5,$6,$7,$8}' | sort -k3 -rn   | grep -v Super-Scaffold_2765
B Super-Scaffold_128 286304 289383 scaffold391 272416 275581 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
B Super-Scaffold_128 284594 284855 scaffold391 270488 270746 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
B Super-Scaffold_128 266292 266691 scaffold391 225960 226359 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
B Super-Scaffold_128 255581 255824 scaffold391 222259 222496 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
B Super-Scaffold_128 250243 253129 scaffold391 217434 221099 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
B Super-Scaffold_128 235682 240996 scaffold391 206111 211508 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 228545 233030 scaffold391 199473 204307 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 196637 199077 scaffold391 166537 168280 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 194367 195477 scaffold391 160376 165420 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 184368 185848 scaffold391 154408 155734 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 180919 181450 scaffold391 149270 149801 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 157395 163162 scaffold391 128066 133869 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 147054 156417 scaffold391 121069 125256 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
B Super-Scaffold_128 125473 131689 scaffold391 99472 105676 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 122495 123949 scaffold391 96401 98153 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 118485 118950 scaffold391 92374 92839 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 97986 100920 scaffold391 71796 74730 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 94753 96023 scaffold391 68568 69838 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 90926 93668 scaffold391 64742 67483 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 86365 87366 scaffold391 56419 57417 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 65228 69299 scaffold391 34479 38677 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 58331 60911 scaffold391 27649 30222 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 46219 51183 scaffold391 15164 20254 Super-Scaffold_128 36603 289383 scaffold391 3965 275581
A Super-Scaffold_128 42639 44843 scaffold391 11559 14161 Super-Scaffold_128 36603 289383 scaffold391 3965 275581

```
See this handdrawn approximation. 

![hand drawn][hand_drawn]

So that means we may need a more sensistive approach for determinig points of recombination than working at the region
level since it could happen within the region level. A=green, B=Orange, p=*E. indica* 

![117 A=green, B=Orange, p=*E. indica* ][117]


[ideal]:ideal_overlap.png
[complex]:complex.png
[complex_sub1]:highlight_acomplex.png
[hand_drawn]:SS128.png
[117]:117_len.png