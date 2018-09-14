"""
This set of functions will link all hits together into a network.
"""

import dc_parse
import dc_links_to_net
import networkx as nx
import matplotlib.pyplot as plt

def make_bed_from_node_name(node_name):
    """

    :param node_name: takes standard node that is bed3 format and returns bed string.
    :return:
    """
    node_name = node_name.split("_")
    chrom = "_".join(node_name[0:-2])
    start = node_name[-2]
    stop = node_name[-1]
    return [chrom,start,stop]

def test_for_over_lap(u, pybed):
    """

    :param u: this is a node name that is back converted into a bed.
    :param pybed: BedTools object that we check for overlap against.
                  this is a bed4 object so it must be explicitly parsed.
    :return: [overlap - boolean ,ret_ent- the single entry that overlaps with u in pybed ]
    """
    overlap = True
    ret_ent = ''
    chrom,start,stop = make_bed_from_node_name(u)
    new_pybed = dc_parse.create_pybed3(chrom,start,stop)
    overlap_test = new_pybed.intersect(pybed, wo=True)
    if len(overlap_test) == 1 :
        overlap = False
        ret_ent = str(overlap_test[0]).split("\t")[3:7]
    return [overlap,ret_ent]

def convert_bed4_raw_ent_to_dict(raw_ent):
    """

    :param raw_ent: moves all bed4 objects into a dictionary for downstream comparison.
                    bed4 objects must be formatted as they are created in dc_parse.py
    :return: {gene:[percent_id,...]}
    """


    ret_dict = {}
    # first strip final _id from
    pers_genes = raw_ent[3] # this should be the fourth
    pers_genes = pers_genes.split("_id")
    assert len(pers_genes) >= 2, "Malformed bed4 comment. It must end in \"id\""
    if len(pers_genes) == 2 : # there was only _id in whole string.
        pers_genes = pers_genes[0]
    elif len(pers_genes) > 2 : #_id occured somewhere in the name.
        pers_genes = "_id".join(pers_genes[0:-1])
    # now split bed4 comment.
    pers_genes = pers_genes.split("^")
    #print (pers_genes[1],pers_genes[2]) # where we have put percent id and genes.
    mean = pers_genes[0]
    ret_dict['mean'] = float(mean)
    pers = [float(x) for x in  pers_genes[1].split(",")]
    genes = [x for x in pers_genes[2].split(",")]
    for per,gene in zip(pers,genes):
        if gene not in  ret_dict :
            ret_dict[gene] = [per]
        else :
            ret_dict[gene].append(per)

    return ret_dict

def gene_wise_compare(n_bed, u_bed, n_dict,u_dict):
    """

    :param n_bed: bed coords in which genes are compared
    :param u_bed:  bed coords in which genes are compared
    :param n_dict:  dict with gene id as key value and [percent id,...] as return value.
    :param u_dict:  same as above.
    :return:  2formatted lists for writing to file per gene and per syn block.
    """

    out_write_per_gene = []
    out_write_by_bed = []
    for gid in n_dict:
        if gid in u_dict:
            if len(n_dict[gid]) == 1 and len(u_dict[gid]) == 1 :
                if n_dict[gid] > u_dict[gid]:
                    out_write_per_gene.append([gid,n_bed,n_dict[gid], u_bed,u_dict[gid]])
                    out_write_by_bed.append([n_bed,'A',gid])
                    out_write_by_bed.append([u_bed,'B',gid])
                elif n_dict[gid] < u_dict[gid]:
                    out_write_per_gene.append([gid,u_bed,u_dict[gid], n_bed,n_dict[gid]])
                    out_write_by_bed.append([u_bed,'A',gid])
                    out_write_by_bed.append([n_bed,'B',gid])
                else: # in case values are equal.
                    out_write_per_gene.append(["equal_"+str(gid),n_bed,n_dict[gid], u_bed,u_dict[gid]])
    return [out_write_per_gene, out_write_by_bed]

def bed_region_compare(u,n,pybed):
    out_write_by_bed = 'Null'
    out_write_per_gene = 'Null'
    if u != n:  # testing for weirdness.
        # now test for overlap.
        u_test, u_bed4 = test_for_over_lap(u, pybed)
        n_test, n_bed4 = test_for_over_lap(n, pybed)
        if u_test is False and n_test is False:
            n_dict = convert_bed4_raw_ent_to_dict(n_bed4)
            u_dict = convert_bed4_raw_ent_to_dict(u_bed4)

            out_write_per_gene, out_write_by_bed = gene_wise_compare(n, u, n_dict, u_dict)
            #for X in out_write_by_bed:
               # print(X)
    return [out_write_per_gene, out_write_by_bed]

def parse_self_vs_parent_for_linkage(tmp_dir,in_file, out_file_prefix):
    """
    first we parse self vs out. This gives us connections between syntenic blocks, based on the
    the out genome. For this to be effective. The syntenic block cutoff must be stringent enough to
    exclude any relict WGDs, this is a particular concern for plants.

    :param tmp_dir:
    :param in_file:
    :param out_file_prefix:
    :return:
    """

    pybed, counter, added_coords = dc_parse.parse_parent_vs_self_syn(stem=in_file)
    """
    Now save pybed to disk as a checkpoint.
    Also save network connection to disk as checkpoint. 
    """
    pybed.saveas('{tmp_dir}/{out_file}.bed'.format(tmp_dir=tmp_dir,out_file=out_file_prefix))
    f = open('{tmp_dir}/{out_file}.net'.format(tmp_dir=tmp_dir,out_file=out_file_prefix), 'w')
    for X, Y in zip(added_coords[0][0::2], added_coords[0][1::2]):
        f.write('{X}\t{Y}\n'.format(X=X, Y=Y))
    f.close()

    """
    Now make a dictionary for all connections. 
    We will want to ref this from b (as parent) to a (as self)
    """
    id_dict = {}
    for b, a in zip(added_coords[0][0::2], added_coords[0][1::2]):
        if b in id_dict:
            id_dict[b].append(a)
        else:
            id_dict[b] = [a]
        #print(b)
    return [id_dict, pybed]


def compare_genes(id_dict, pybed):
     final_out_per_gene = []
     final_out_per_bed = []
     tested = []
     linker = 'Null'
     for b_region_id in id_dict:
         # print(b_region)
         out_write_per_gene, out_write_by_bed = ['Null', 'Null']
         if b_region_id not in tested:
             print(b_region_id) # all 3 are correctly spit out.
             chrom, start, stop = make_bed_from_node_name(b_region_id)
             b_region = dc_parse.create_pybed3(chrom, start, stop)
             overlaps = b_region.intersect(pybed, wo=True, )
             print(len(overlaps))
             tested.append(b_region_id)

             # this block limits number of hits to 2.
             # it uses default bedtools intersect so it exlcudes any
             # overlaps. This helps avoid ambiguity because we are looking
             # for solid anchor points.

             if len(overlaps) == 1:
                 linker = b_region_id
             elif len(overlaps) == 2:
                 linker = "_".join(str(overlaps[0]).split("\t")[3:6])
             else:
                 linker = 'Null'  # here we are enforcing no overlaps among regions. If this is too strict, we can start
                 # start adding some more subtle parsing here.
             # get self region
             # check self region against bed if no other overlaps use.
             if linker != 'Null':
                 if linker == b_region_id:
                     # this is the case where bed coordinates perfectly overlap.
                     if len(id_dict[b_region_id]) == 2:
                         u, n = id_dict[b_region_id]
                         out_write_per_gene, out_write_by_bed = bed_region_compare(u, n, pybed)
                 elif linker != b_region_id:
                     if len(id_dict[linker]) == 1 and len(id_dict[b_region_id]) == 1:
                         u = id_dict[b_region_id][0]
                         n = id_dict[linker][0]
                         out_write_per_gene, out_write_by_bed = bed_region_compare(u, n, pybed)
         if out_write_per_gene != 'Null' and out_write_by_bed != 'Null':
             final_out_per_gene += out_write_per_gene
             final_out_per_bed += out_write_by_bed
     return [final_out_per_bed, final_out_per_gene]

def write_final_files(final_out_per_bed,final_out_per_gene,tmp_dir,out_file_prefix):
    """

    :param final_out_per_bed: formatted listed returned by compare_genes()
    :param final_out_per_gene: ''               ''
    :param tmp_dir: directory to write file to.
    :param out_file_prefix: unique tag for output.
    :return:
    """
    f_bed_out = open('{tmp_dir}/syn_block_ab_gene_{out_file_prefix}.txt'.format(tmp_dir=tmp_dir,
                                                                                 out_file_prefix=out_file_prefix)
                     , 'w')
    for X in final_out_per_bed:
        f_bed_out.write("\t".join(X))
        f_bed_out.write("\n")
    f_bed_out.close()

    f_gene_out = open('{tmp_dir}/gene_asyn_aper_bsyn_bper_{out_file_prefix}.txt'.format(tmp_dir=tmp_dir,
                                                                                         out_file_prefix=out_file_prefix
                                                                                         ), 'w')
    for X in final_out_per_gene:
        f_gene_out.write("\t".join([str(Y) for Y in X]))
        f_gene_out.write("\n")
    f_gene_out.close()

def ab_syn_finder_main(tmp_dir,out_file_prefix, in_file):
    """

    :param tmp_dir: this is the directory you want analysis written to. Default will be wd.
    :param out_file_prefix: common string to include in the name of all outfiles.
    :param in_file: DAG_Chainer output.
    :return: Writes to disk preliminary network that shows basic network recovered from simple parsing of dc output.
             It contains no connections discovered through pybedtools.intersect (.net)
             writes bed of of all syntenic regions to be parsed from dc output (.bed).
             syn_block_ab_gene{out_file_prefix}.txt contains syntenic block, AvsB call, gene
             gene_asyn_aper_bsyn_bper{out_file_prefix}.txt
             gene, predicted A syntenic block, percent identity shared with parent,\
             predicted B syntenic block, percent identity shared with parent
             Note: if equal_ precedes gene name percent ids were identical.
    """

    id_dict, pybed = parse_self_vs_parent_for_linkage(tmp_dir,in_file,out_file_prefix)

    pybed3 = dc_links_to_net.pybed4_convert_pybed3(pybed)
    intersect_bed = pybed3.intersect(pybed3,
                                     wo=True)
    syn_graph = nx.Graph()
    for idp in id_dict:
        for ida in id_dict[idp]:
            if syn_graph.has_edge(idp,ida) is False:
                syn_graph.add_edge(idp,ida)
                syn_graph.node[idp]['name'] = 'p' # for parent genome
                syn_graph.node[ida]['name'] = 'a' # for the ab genome candiate.

                # this is an undirected graph.

    syn_graph = dc_parse.add_intersect_to_syn_graph(syn_graph, intersect_bed)
    labels = dict((n,d['name']) for n,d in syn_graph.nodes(data=True))
    nx.draw(syn_graph,labels=labels)
    plt.show()
    for sg in nx.connected_component_subgraphs(syn_graph):
        print(type(sg))
        cyc = nx.cycle_basis(sg.to_undirected())
        print(cyc)
        if len(cyc) ==0 :
            parent_nodes =[]
            for n in sg.node :
                if sg.node[n]['name'] == 'p':
                    print("found_nodes",sg[n],n)
                    print(type(n))
                else :
                    pass


                print(n)

            #eb = nx.edge_boundary(sg,)
            #print(eb)

            bed_list = []
            print(sg)
            for node in sg:
                bed_list.append(make_bed_from_node_name(str(node)))
            print(bed_list)
    #final_out_per_bed, final_out_per_gene =  compare_genes(id_dict=id_dict, pybed=pybed)
    #write_final_files(final_out_per_bed,final_out_per_gene,tmp_dir,out_file_prefix)


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
    argv = parser.parse_args()

    ab_syn_finder_main(argv.out_dir,argv.out_pre, argv.in_file)












