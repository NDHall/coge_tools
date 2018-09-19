##CoGe steps
1. select parent and child genomes.
2. generate ![synmap][coge_syn1]
    Here we are using ks/kn calculations in tandem with quota align 1(*E. coracana*) : 2(*E. indica*). Initially, 
    I chose to avoid quota align(qac) because when using a gene based criteria, it is possible to get physical 
    overlaps on the chromosomes. These overlaps sometimes lead to double dipping during gene analysis, such that one gene
    may be used for 2 different syntenic blocks. Though, my work with the qac files suggests this happens, it seems to 
    be rare enough that we will be able to simply filter out any gene that appears more than twice. If you would really
    like to avoid a 'physical' overlap this can be done, in theory, by tweaking qac and SynMap parameters. For example 
    choosing to use coordinates instead of cds to determine syntenic blocks would be a practical first step. 
3. check SynMap qac output. One thing that we have noticed is that even though a ks value may be computed for all
    positions it is reported for a subset of positions. Also, ks can become very high. This is because the underlying
    codeml program used to call ks does not put an upward limit on ks values. Though in practice we are looking for ks 
    values that are relatively low since we expect the divergence between *E. indica*, *E. coracana* A and B to be
    small and recent. We have been able to show this in past *E. coracana* vs self ks comparisons. That show divergence 
    between A and B genomes to have about a maximum  mean of about 0.5*ks* per syntennic block, wit the signatures of 
    ancient WGDs occurring after this cutoff.
    
##DAGChainer Filtering
We want to now filter the output so that we return genes or anchors from syntenic blocks that pass filtering cutoffs and
only occur twice. Since these are closely related a minimum average percent identity of 90 is required for
each syntenic block. 
    
4. If you have a a good sense of where the cutoff exists in the data, or if you have used qac, calling by percent id 
    alone will yield more hits. Here we compare the use of `strict_ks = False` to `strict_ks = True`.
    ```bash
       #qac == False
       ndh0004@IorekByrnison:~/code/coge_tools/test_out$ awk '{print $4}' bog_sep18_abcalls.tsv |\
                                                         sort | uniq -c |awk '{print $1}' | sort | uniq -c | sort -n 
       10409 2 # number of anchors that occur 2 times
    
       #qac == True
    (venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ awk '{print $4}'  51576_52024_qac_ks_sep19_abcalls.tsv | sort | uniq -c |awk '{print $1}' | sort | uniq -c | sort -n 
       6398 2
 
       #qac == False and strict_ks == True
    ndh0004@IorekByrnison:~/code/coge_tools/test_out$ awk '{print $4}' bog_sep18_strict_abcalls.tsv |\
                                                       sort | uniq -c |awk '{print $1}' | sort | uniq -c | sort -n 
       2515 2 # number of anchors that occur 2 times
      ndh0004@IorekByrnison:~/code/coge_tools/test_out$
   
5. Now run ab_caller.py 





[coge_syn1]:(https://genomevolution.org/coge/SynMap.pl?dsgid1=51576;dsgid2=52024;D=20;A=5;w=0;b=6;ft1=1;ft2=1;autogo=1;Dm=0;tdd=10;gm=0;snsd=0;ma=1;da=1;do1=2;do2=1;do=50;fb_ws=100;fb_nqc=25;fb_ntc=25;fb_rru=1;cs=1;cmin=0;cmax=3;logks=0;sr=1;cso=S;dt=geneorder;ks=3)