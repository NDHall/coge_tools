###Pipeline for using dc_tools.

[CoGe SynMap][dc_link] produces data rich output in the form of DAGChainer
files. Scripts from **dc_tools** can be chained together to 
analyze the output and search for A vs B regions of genome. The 
goal in this case is to determine the data are doing. 

####Data Visualization

1. Begin with raw [DAGChainer][dc_link] (these are pretty far down under results) files. 
2. **run dc_gene_count.py** Create histogram of genes per syntenic block distribution, and use this to set your cutoffs.
 For example, if you 
examined diploid parent vs diploid parent syntenic block, it will show syntenic regions created during a WGD. These 
blocks would be the minimum block size when looking for A vs B regions, so that all known WGD blocks are screened by 
default.
3. **run dc_links_to_net.py** This will show how syntenic regions overlap. Remember, SynMap is looking at gene order 
between genomes. This can lead to overlapping syntenic blocks that will create ambiguity for us. Using this script will 
help determine how your results are overlapping. I suggest using a network visualization tool like [gephi][gephi], or digging 
into the sub-graphs of interest with your own custom method. [Here][example1] is example code I used in my approach.
####ReRun SynMap 
4. Adjust SynMap parameters. Now that we know if regions are overlapping, we can look critically at syntenic blocks and 
set SynMap parameters so that it calls more conservatively, than its default behavior. 
####Downstream Analysis
5. **run ab_syn_finder.py** on new DAGChainer output. 
6. **run ab_call.py** using output from **ab_syn_finder.py** to make calls for each syntenic region.
7. Examine output and look at ambiguous calls. Ambiguous calls will either be for a series of low quality regions or, if 
you're lucky, a homelogous cross over event. Your confidence in this assesment will be based on your confidence in your 
over all assemblies, and the paramenters you have set before. 




[dc_link]:https://genomevolution.org/wiki/index.php/SynMap
[gephi]:https://gephi.org/
[example1]:https://github.com/NDHall/coge_tools/blob/master/dc_tools/test_edges.py
