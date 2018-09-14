import networkx as nx
import matplotlib.pyplot as plt
import ab_syn_finder
import dc_parse
import dc_links_to_net
import pybedtools


def create_gene_dict(pybed4,id_dict):
    ret_dict = {}
    for entry in pybed4:
        chrom = entry[0]
        start = entry[1]
        stop = entry[2]
        id = '{chrom}_{start}_{stop}'.format(
            chrom=chrom,
            start=start,
            stop=stop)
        if id not in id_dict: # first only take A genome values.
            if id not in ret_dict: # now we only take unique values.
                ret_dict[id] = ab_syn_finder.convert_bed4_raw_ent_to_dict(entry) # now we have genes linked to A genome
                #ids
            else:
                # in this case we have repeated region in the child genome. We don't want it because it is
                # now a degenerate comparison.
                ret_dict[id] = None
    return ret_dict


def get_connected_target_blocks(syn_graph):
    ret_anodes = []

    for sg in nx.connected_component_subgraphs(syn_graph):
        cyc = nx.cycle_basis(sg.to_undirected())
        anodes = []
        if len(cyc) ==0 :
            clean_graph = True

            for n in sg.node :
                if n not in anodes and sg.node[n]['name'] == 'a':
                    anodes.append(str(n))
                source_list =[sg.node[n]['name']]
                for neighbor in sg.neighbors(n):
                   source_list.append(sg.node[neighbor]['name'])
                if 'a' not in source_list:
                    #no set of all p nodes is acceptable.
                    clean_graph = False
                if source_list.count('a') > 2 :
                    # 2 is okay because we can have cases where 1 p node links 2 a nodes.
                    clean_graph =False
                if source_list[0] == 'a' \
                    and len(source_list) != 2:
                     # all nodes that start a must only have 1 edge.
                    clean_graph = False
                if 'p' not in source_list:
                    # all a nodes must link only through p.
                    clean_graph = False
            if clean_graph == True:
                ret_anodes.append(anodes)
    return ret_anodes


def bag_of_genes(gene_by_region_dict, list_of_grouped_nodes):
    genes_to_call = {}
    for group in list_of_grouped_nodes:
        if len(group) >1: # only take good nodes.
            for node in group:
                for gene in gene_by_region_dict[node]:
                    if gene != 'mean' and gene is not None:  # we don't need to add all the means to the dict
                        if len(gene_by_region_dict[node][gene]) == 1: # make sure no weirdness happened in parsing.
                            per_id = gene_by_region_dict[node][gene][0]
                            if gene not in genes_to_call:
                                genes_to_call[gene] = [[node,per_id]]
                            else:
                                genes_to_call[gene].append([node,per_id])
    return genes_to_call





def call_genes(out,bag_of_genes):
    fout = open('{out}_abcalls.tsv'.format(
        out=out
    ), 'w')
    fout_no = open('{out}_nocalls.tsv'.format(
        out=out
    ), 'w')
    for gene in bag_of_genes:
        if len(bag_of_genes[gene])== 2:
            n,npid = bag_of_genes[gene][0]
            u, upid = bag_of_genes[gene][1]
            if npid > upid :
                ncall = 'A'
                ucall = 'B'
            elif upid > npid:
                ncall = 'A'
                ucall = 'B'
            else:
                ncall = 'eq'
                ucall = 'eq'
            for call,region,pid in [[ncall,n,npid],[ucall,u,upid]]:
                sout = '{r}\t{p}\t{c}\t{g}\n'.format(
                    r=region,
                    c=call,
                    p=pid,
                    g=gene
                )
                fout_no.write(sout)
            if ncall != 'eq':
                nout_string = '{block}\t{call}\t{gene}\n'.format(
                    block=n,
                    call=ncall,
                    gene=gene
                )
                uout_string= '{block}\t{call}\t{gene}\n'.format(
                    block=u,
                    call=ucall,
                    gene=gene
                )
                for s in [nout_string,uout_string]:
                    fout.write(s)
        else:
            for region,pid in bag_of_genes[gene]:
                sout = '{r}\t{p}\tn{c}g\t{g}\n'.format(
                        r=region,
                        c=len(bag_of_genes[gene]),
                        p=pid,
                        g=gene
                    )
                fout_no.write(sout)
    fout_no.close()
    fout.close()



def main_bag_of_genes(out_dir,out_file,infile, mean_cutoff):
    id_dict, pybed4 = ab_syn_finder.parse_self_vs_parent_for_linkage(out_dir, infile, out_file)
    # parse out gene by region info.
    gene_by_region_dict = create_gene_dict(pybed4,id_dict)
    pybed3 = dc_links_to_net.pybed4_convert_pybed3(pybed4)
    intersect_bed = pybed3.intersect(pybed3,
                                     wo=True)
    syn_graph = nx.Graph()
    for idp in id_dict:
        for ida in id_dict[idp]:
            mean = gene_by_region_dict[ida]['mean'] # we should not not need to check if this is in dict.
            if syn_graph.has_edge(idp,ida) is False and \
                    mean >= mean_cutoff: # here we can maximize the number of
                                         # comparisons by throwing out any connections
                                         # that fall below a certain quality.
                syn_graph.add_edge(idp,ida)
                syn_graph.node[idp]['name'] = 'p' # for parent genome
                syn_graph.node[ida]['name'] = 'a' # for the ab genome candiate.
                # this is an undirected graph.
    syn_graph = dc_parse.add_intersect_to_syn_graph(syn_graph, intersect_bed)
    # now lets get the nodes that are connected by strong graph shapes.
    connected_anodes =  get_connected_target_blocks(syn_graph)
    # now for each set of connected nodes lets use a bag of genes approach.
    bogenes = bag_of_genes(gene_by_region_dict,connected_anodes)
    call_genes(out_file,bogenes )



def message(case):
    if case == 'help':
        string_ret = """
        ab_syn_finder.py takes DAG_Chainer(dc) output from CoGe in which a tetraploid 
        genome has been compared to a known diploid parent. The input file should 
        be organized such that tetraploid parent is in the a column of the dc output
        and the diploid parent is in the b column.  
        """

if __name__ == '__main__':
    import argparse
    import os
    cwd = str(os.getcwd())

    parser = argparse.ArgumentParser(description=message('help'))
    parser.add_argument('-i','--infile',
                        help='dc file, tetraploid genome in a column and diploid genome in b column',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-d','--dir',
                        help='directory to save results to, default is current working directory',
                        required=False,
                        default=cwd,
                        type=str,
                        dest='out_dir')
    parser.add_argument('-o','--outfile',
                        required=True,
                        type=str,
                        dest="out_pre")
    parser.add_argument('-m','--mean_perid_cutoff',
                        help='minimum mean percent id to consider default==90%',
                        required=False,
                        default=90.0,
                        type=float,
                        dest='mean_cutoff')
    argv = parser.parse_args()

    main_bag_of_genes(argv.out_dir,argv.out_pre, argv.in_file, argv.mean_cutoff)


