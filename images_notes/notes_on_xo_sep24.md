###Notes on Cross over Chrom Calls.
I have some concerns about calling  cross over regions from the syntentic region call. 
These will extend to general AB calls as well, but by the end of these notes, I hope to provide some 
insight into the best filtering practices. The main problem is largely there are AB cross overs that can be 
dropped because of ambig call, but likely represent a true crossing over.

```bash
#list of AB calls provisional and strongndh0004@IorekByrnison:~/code/coge_tools/test_out$ egrep '^AB'  bog_sep24_qaccalled_chrom_calls_all.tsv 
AB	provisional	14	3	Super-Scaffold_50	220	65	9	15
AB	provisional	7	1	Super-Scaffold_25	126	27	4	1
AB	provisional	12	1	Super-Scaffold_70	104	115	26	20
AB	strong	2	0	scaffold591_size730508_subseq_1:450970_obj	11	26	0	0
AB	strong	4	0	Super-Scaffold_13	41	43	0	0
AB	provisional	20	4	Super-Scaffold_84	225	250	30	37
AB	provisional	5	2	Super-Scaffold_291	29	34	7	3
AB	provisional	3	1	scaffold306_size1284463_obj	38	44	5	5
AB	provisional	5	1	Super-Scaffold_210	65	13	4	11
AB	strong	4	0	Super-Scaffold_105	21	7	0	0
AB	provisional	12	2	Super-Scaffold_95	187	24	20	12
AB	strong	12	0	Super-Scaffold_2412	35	225	0	0
AB	strong	6	0	Super-Scaffold_146	19	72	0	0
AB	provisional	20	1	Super-Scaffold_39	65	292	12	19
AB	provisional	6	1	Super-Scaffold_266	24	71	1	5
AB	strong	9	0	Super-Scaffold_311	134	31	0	0
AB	provisional	5	1	Super-Scaffold_76	37	13	1	0
AB	strong	7	0	Super-Scaffold_307	25	83	0	0
AB	strong	2	0	scaffold93_size1230182_subseq_1:987274_obj	14	15	0	0
AB	strong	14	0	Super-Scaffold_12	65	220	0	0
AB	provisional	2	1	Super-Scaffold_2765	36	18	5	6
AB	provisional	13	1	Super-Scaffold_10	114	126	0	1
AB	strong	12	0	Super-Scaffold_128	27	297	0	0
AB	strong	10	0	Super-Scaffold_253	150	20	0	0
AB	strong	4	0	Super-Scaffold_1173	43	41	0	0
AB	provisional	3	2	Super-Scaffold_202	5	17	9	2
AB	provisional	6	2	Super-Scaffold_158	71	49	26	19
AB	provisional	4	1	Super-Scaffold_147	12	50	11	4


```

##Cross over site provisional example. 
![Super-Scaffold_10 shown in igv. B calls == blue, A calls == purple, gray == all][ss_50_xo]

Let's look at the linked scaffold that is the the likely site of crossing over. 

![Super-Scaffold_1693 shown in igv. B calls == blue, A calls == purple][ss_1693]



Now cross over site and upstream B, on SS_50 (not shown) are both called ambig, but represent a 
a real possibility for crossing over. 
```bash

ndh0004@IorekByrnison:~/code/coge_tools/test_out$ grep Super-Scaffold_1693 bog_sep24_qaccalled_region_calls.tsv 
ambig	Super-Scaffold_1693	2918	135878	9	7	0.6170750774519739	0.25

ndh0004@IorekByrnison:~/code/coge_tools/test_out$ egrep -w Super-Scaffold_50 bog_sep24_qaccalled_region_calls.tsv| awk '$4<1200000'  
ambig	Super-Scaffold_50	944048	1128041	7	9	0.6170750774519739	0.25
ambig	Super-Scaffold_50	40601	96331	2	5	0.2568392579578533	1.2857142857142858
ndh0004@IorekByrnison:~/code/coge_tools/test_out$ 

#566 is the scaffold that is the other call for ss_50

ndh0004@IorekByrnison:~/code/coge_tools/test_out$ egrep -w Super-Scaffold_566 bog_sep24_qaccalled_region_calls.tsv  A	Super-Scaffold_566	852770	902763	7	1	0.033894853524689295	4.5
ambig	Super-Scaffold_566	1339186	1505534	5	2	0.2568392579578533	1.2857142857142858

```

###Cross over site strong 

Not only is it reciprocal, it is through only 1 *Eleusine indica* contig.
![Connection scheme][ss_13_xo_con]

![Super-Scaffold_13 shown in igv. B calls == blue, A calls == purple][ss_13_xo]
It is worth noting that recip. calls will not always end up so nicely aligned between scaffolds.
![Super-Scaffold_1173 shown in igv. B calls == blue, A calls == purple][ss_1173_xo]


###Interesting possible cross over

![Super-Scaffold_39 shown in igv. B calls == blue, A calls == purple][ss_39_xo]
![Super-Scaffold_95 shown in igv. B calls == blue, A calls == purple][ss_95_xo]

It does not look like this the result of mis-assembly since we see no breaks in contigs.








[ss_1693]:ss_1693.png
[ss_50_xo]:ss_50_xo.png
[ss_13_xo_con]:ss_13_xo_con.png
[ss_13_xo]:ss_13_xo.png
[ss_1173_xo]:ss_1173_xo.png
[ss_39_xo]:ss_39_xo.png
[ss_95_xo]:ss_95_xo.png