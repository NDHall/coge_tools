"""
dc provides links between 2 syntenic blocks.
ab_syn_finder outputs list of links. The purpose.
is to use either list output from ab_syn_finder.
or a raw dc to produce a network representation of
how syntenic regions are connected.
ToDo: update to take previous lists. currently it only takes raw dc output
"""

import dc_parse
import dc_gene_count
import networkx as nx
import pybedtools
import matplotlib.pyplot as plt

def create_and_save_net_and_bed(in_file,out_file_prefix,tmp_dir):
    pybed, counter, added_coords = dc_parse.parse_parent_vs_self_syn(stem=in_file)
    """
    Now save pybed to disk as a checkpoint.
    Also save network connection to disk as checkpoint. 
    """
    pybed.saveas('{tmp_dir}/{out_file}.bed'.format(tmp_dir=tmp_dir,out_file=out_file_prefix))
    f = open('{tmp_dir}/{out_file}.net'.format(tmp_dir=tmp_dir,out_file=out_file_prefix), 'w')
    syn_graph=nx.Graph()
    for X, Y in zip(added_coords[0][0::2], added_coords[0][1::2]):
        f.write('{X}\t{Y}\n'.format(X=X, Y=Y))
        if syn_graph.has_edge(X,Y) is False:
            syn_graph.add_edge(X,Y)

    f.close()
    #print(len(syn_graph))
    return [pybed,syn_graph]

def pybed4_convert_pybed3(pybed4):
    pybed_3 = "\n".join(['\t'.join(str(bed).split("\t")[0:3]) for bed in pybed4 ])
    pybed_return = pybedtools.BedTool(pybed_3,from_string=True)
    return pybed_return

def main_links_to_net(show,dpi,bins,tmp_dir, out_file_prefix,in_file,graphic_type):

    #in_file = '/home/ndh0004/code/coge_tools/test_data/SScaf_indica_v_cor.test'
    #out_file_prefix = 'CvI_ss10'
    #graphic_type = 'png'

    pybed, syn_graph = create_and_save_net_and_bed(in_file,out_file_prefix,tmp_dir)
    pybed3 = pybed4_convert_pybed3(pybed)
    #print(pybed3)
    intersect_bed = pybed3.intersect(pybed3,
                                     wo=True)
    syn_graph = dc_parse.add_intersect_to_syn_graph(syn_graph,intersect_bed)
    #print(len(syn_graph))
    nx.draw(syn_graph)
    plt.savefig('{tmp_dir}/{out_file}.{graph_type}'.format(tmp_dir=tmp_dir,out_file=out_file_prefix,
                                                                 graph_type=graphic_type))
    nx.write_edgelist(syn_graph,'{tmp_dir}/{out_file}.{graph_type}'.format(tmp_dir=tmp_dir,out_file=out_file_prefix,
                                                                 graph_type='edges'))

    if show is True:
        plt.show()

    sub_graph_len = [len(con_el) for con_el in nx.connected_components(syn_graph)]


    out = '{tmp_dir}/{out_file}.{graph_type}'.format(tmp_dir=tmp_dir, out_file=out_file_prefix,
                                                     graph_type='hist.' + graphic_type)
    dc_gene_count.make_histogram(sub_graph_len, bins,out,show,dpi )


if __name__ == '__main__':
    # parse raw_dc
    import argparse
    import os

    cwd = str(os.getcwd())

    parser = argparse.ArgumentParser(description="""
    Provides a network of showing connections among syntenic blocks and a histogram indicating frequency of sub-graph
    length occurences. If it is and A vs AB genome 3 would be an ideal sub-graph length.
    """)

    parser.add_argument('-i', '--infile',
                        help='dc file',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-d', '--dir',
                        help='directory to save results to, default is current working directory',
                        required=False,
                        default=cwd,
                        type=str,
                        dest='out_dir')
    parser.add_argument('--dpi',
                        help='dpi for output graph',
                        required=False,
                        default=300,
                        type=int,
                        dest='dpi')
    parser.add_argument('-g','--graphic_type',
                        help='extension indicating type of graphic output default png',
                        required=False,
                        default='png',
                        type=str,
                        dest='graphic_type')
    parser.add_argument('--bins',
                        help='size of bins for histogram',
                        required=False,
                        default=20,
                        type=int,
                        dest='bins')
    parser.add_argument('-s', '--show',
                        help='show graph when made',
                        required=False,
                        default=False,
                        dest='show',
                        action="store_true")
    parser.add_argument('-o', '--outfile',
                        required=True,
                        type=str,
                        dest="out")

    argv = parser.parse_args()

    out = '{dir}/{out}'.format(dir=argv.out_dir, out=argv.out.rstrip("/"))

    main_links_to_net(show=argv.show,
                      bins=argv.bins,
                      dpi=argv.dpi,
                      tmp_dir=argv.out_dir,
                      out_file_prefix=argv.out,
                      in_file=argv.in_file,
                      graphic_type=argv.graphic_type
                      )





