###Notes on ab_calls.py
After successful run on on hopper. I found that we had uneven numbers of scaffolds called 
for provisional. This is not unexpected given that scaffolds are not the unit from which region
calls are made. However, it is worth confirming that all has behaved as expected. 

<provisional_calls>

|call |call_qual  |num_good_reg|num_ambig_reg|chrom                 |A_gene_calls|B_gene_calls|ambig_a|ambig_b|
|:---:|:---------:|:------------|-------------|--------------------- |------------|------------|-------|-------|
|AB   |provisional|8           |1            |Super-Scaffold_50     |145         |54           |7     |9      |
|AB|provisional|3|1|Super-Scaffold_210|37|12|4|11|
|AB|provisional|3|1|Super-Scaffold_147|12|37|11|4|
|AB|provisional|3|1|Super-Scaffold_158|44|16|18|17|
|B|provisional|1|1|Super-Scaffold_266|0|18|9|10|
|A|provisional|1|1|Super-Scaffold_100|18|0|10|9|
|AB|provisional|6|1|Super-Scaffold_39|37|74|12|19|
|AB|provisional|4|1|Super-Scaffold_95|56|20|19|12|
|B|provisional|8|2|**Super-Scaffold_84**|74|115|28|28|

notice that Super-Scaffold_84 has 2 ambig calls.  We weant to verify that these are both to scaffolds already in 
provisional. This is because you cannot have less than 2 ambig calls. If one syntenic region is ambig its mate by 
defintion is also ambig. But since scaffolds are different sizes it means that we have to also consider the ambig calls
in tandem with provisional calls consider this situation.

![alt text][ambig_syn]

In this case we would get 1 provisional, 1 strong and 1 ambig at the chrom level but at the region level we simply have
2 ambig calls and 2 strong calls. 

We can verify that our calls are produced in 2s and that our network is correct by searching the chrom call files. 
```bash
(venv) [ndh0004@node019 full_run_sept5]$ grep -A1 -B1 ambig gt4_ab_calls_region_calls.tsv 
B	Super-Scaffold_243	1470117	1573006	0	16	6.334248366623988e-05	16.0
ambig	Super-Scaffold_147	2303797	2424319	11	4	0.07070114486598289	3.2666666666666666
ambig	Super-Scaffold_210	2719518	2936984	4	11	0.07070114486598289	3.2666666666666666
A	Super-Scaffold_214	600747	861824	15	3	0.004677734981047276	8.0
--
B	Super-Scaffold_2591	142640	364699	1	24	4.224909405005696e-06	21.16
ambig	Super-Scaffold_266	2601345	2877809	9	10	0.8185458083820435	0.05263157894736842
ambig	Super-Scaffold_100	3448817	3750563	10	9	0.8185458083820435	0.05263157894736842
A	Super-Scaffold_100	3907643	4064447	18	0	2.2090496998585475e-05	18.0
--
B	Super-Scaffold_291	8426441	8580486	0	19	1.3071845366762988e-05	19.0
ambig	Super-Scaffold_158	1181688	1508814	18	17	0.8657723749926214	0.02857142857142857
ambig	Super-Scaffold_302	489702	755453	17	18	0.8657723749926214	0.02857142857142857
A	Super-Scaffold_311	1635514	1858554	20	1	3.381272563301577e-05	17.19047619047619
--
B	Super-Scaffold_176	1821987	2072837	1	17	0.000162440845186801	14.222222222222221
ambig	Super-Scaffold_39	5753501	6115948	12	19	0.2086677876982606	1.5806451612903225
ambig	Super-Scaffold_95	1896656	2247103	19	12	0.2086677876982606	1.5806451612903225
A	Super-Scaffold_95	2384454	2613603	18	1	9.616588268883469e-05	15.210526315789474
--
A	Super-Scaffold_5	935007	1062886	14	3	0.007632881787792295	7.117647058823529
ambig	Super-Scaffold_1693	2918	135878	9	7	0.6170750774519739	0.25
ambig	Super-Scaffold_50	944048	1128041	7	9	0.6170750774519739	0.25
A	Super-Scaffold_50	2256035	2447835	15	0	0.00010751117672950066	15.0
--
B	Super-Scaffold_692	2274	209177	1	13	0.0013406411172294807	10.285714285714286
ambig	Super-Scaffold_7	9892374	10061294	1	0	0.31731050786291115	1.0 # this is the only ss7 call and  so 
ambig	Super-Scaffold_7	3720209	3889129	0	1	0.31731050786291115	1.0    # is discarded
A	Super-Scaffold_70	5747388	6069511	35	0	3.2970532689972886e-09	35.0
--
B	Super-Scaffold_84	3509449	3767784	1	29	3.186355462846913e-07	26.133333333333333
ambig	Super-Scaffold_84	2602492	2884299	9	19	0.05878172135535891	3.5714285714285716 # target scaffold produced 
ambig	Super-Scaffold_84	4040134	4273865	19	9	0.05878172135535891	3.5714285714285716 # both ambigs
A	Super-Scaffold_84	6235900	6484426	22	1	1.193331081541827e-05	19.17391304347826
(venv) [ndh0004@node019 full_run_sept5]$ logout

```




[ambig_syn]:/home/ndh0004/code/coge_tools/images_notes/syn-block.jpg 

