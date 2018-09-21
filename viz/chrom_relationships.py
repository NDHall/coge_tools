import dc_tools.dc_parse
import dc_tools.dc_ks_ka_to_bed as dc_ks_ka_to_bed
import dc_tools.dc_ks_bagofgenes as bog
from data_get import read_bag_of_genes_classy
import networkx as nx
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import Bio.Alphabet as Alphabet


class SpoofSeq(Seq):
    """
    We need to put length into gff. One way would be to define length at the level of the Seq object,
    Which is where WriteGff() is taking the length from at the top of each section. So we need to create
    A Seq object that we can spoof the length on and set it to some number. That makes sense.
    """
    def __init__(self, data, alphabet=Alphabet.generic_alphabet):
        super().__init__( data, alphabet=Alphabet.generic_alphabet)
        self.length = 90000000 # just an arbitrary number to check if it works.

    def __len__(self):
        """

        :return: This is so that we can set and spoof the length of each sequence without having to read the sequence
        into memory. This will be helpful when dealing with whole genomes. As awe are in this case.
        """
        return self.length



def gene_add(classy_gene, syn_graph):
    p_stem = classy_gene.ptregion.split("\t")[0]
    c_stem = classy_gene.ctregion.split("\t")[0]
    pid = 'p_{c}'.format(c=p_stem)
    cid = '{call}_{chrom}'.format(call=classy_gene.call,
                                  chrom=c_stem)
    if syn_graph.has_node(pid) is False:
        syn_graph.add_node(pid,name='p')
    if syn_graph.has_node(cid) is False:
        syn_graph.add_node(cid, name=classy_gene.call)
    if syn_graph.has_edge(pid,cid) is False:
        syn_graph.add_edge(pid,cid,weight=1)
    else:
        syn_graph[pid][cid]['weight'] += 1
    return syn_graph

def chrom_linkages_through_parent(bog,out,filter=True):
    """

    :param bog: takes a bag of genes a list made of class Gene_in_bag_objects:
    :param out string of gephi output file.
    :return: networkx network of chrom linkages for general inspection.
    """

    syn_graph = nx.Graph()

    for classy_gene in bog:
       if filter == True:
           if classy_gene.call == 'A' or \
               classy_gene.call == 'B':
               syn_graph = gene_add(classy_gene,syn_graph)
       else:
           syn_graph = gene_add(classy_gene,syn_graph)


    nx.write_gexf(syn_graph,out)


def chrom_linkages_by_list(bog,out,chrom_list,filter=True,sub_graph=True):
    """

    :param bog: takes a bag of genes a list made of class Gene_in_bag_objects:
    :param out string of gephi output file.
    :return: networkx network of chrom linkages for general inspection.
    """

    syn_graph = nx.Graph()

    for classy_gene in bog:
        if classy_gene.ctregion.split("\t")[0] in chrom_list or \
            classy_gene.ptregion.split("\t")[0] in chrom_list :
            if filter == True:
               if classy_gene.call == 'A' or \
                   classy_gene.call == 'B':
                   syn_graph = gene_add(classy_gene,syn_graph)
            else:
                syn_graph = gene_add(classy_gene,syn_graph)


    nx.write_gexf(syn_graph,out)

    if sub_graph == True :
        counter =1
        for X  in [con for con in nx.connected_components(syn_graph) ]:
            g = syn_graph.subgraph(X)
            nout = out.replace('.gexf' , '_{l}_{c}.gexf'.format(l=len(X),c=counter))
            nx.write_gexf(g,nout)

def create_or_append_gff_feature(classy_gene,gff_dict,chrom_len_dict):
    """

    :param classy_gene: This is class Gene in Bag object
    :param gff_dict: dictionary for creating final gff once everything is processed.
    :param chrom_len_dict: dictionary for determing lenghth of sequence.
    :return:dictionaries that will be used for creating gff file.
    """
    ct_chrom, ct_start, ct_stop = classy_gene.ctregion.split('\t')
    seq = SpoofSeq("ACTG")
    rec = SeqRecord(seq, ct_chrom)
    qualifiers = {"source": "prediction", "score": 10.0, "other": [classy_gene.ptregion.replace('\t', '_')],
                  "ID": classy_gene.ctregion.replace('\t', '_')}

    # Since we are just looking at the region the strand here is set as 1. This could be modified to include the
    # orientation from block, if block gets parsed. However, for these qac.ks files, that we are using, that
    # information does not make it to this stage.

    top_feature = SeqFeature(FeatureLocation(int(ct_start) + 1, int(ct_stop), ),
                             strand=+1,
                             type='region',
                             qualifiers=qualifiers)
    sub_qualifiers = {"source": "prediction", "call": classy_gene.call, "linking_gene": classy_gene.pgene}

    cg_chrom, cg_start, cg_stop = classy_gene.cgene_region.split("\t")
    top_feature.sub_features = [SeqFeature(FeatureLocation(int(cg_start) +1, int(cg_stop)),
                                           type='gene',
                                           strand=classy_gene.cgene_orient[0],
                                           qualifiers=sub_qualifiers)]
    if ct_chrom not in gff_dict:
        gff_dict[ct_chrom] = [rec, top_feature]
    else:
        gff_dict[ct_chrom].append(top_feature)
    if ct_chrom not in chrom_len_dict:
        chrom_len_dict[ct_chrom] = [int(ct_start), int(ct_stop)]
    else:
        chrom_len_dict[ct_chrom] += [int(ct_start), int(ct_stop)]

    return [gff_dict,chrom_len_dict]



def convert_gff(bog,out_file):
    chrom_len_dict = {}
    gff_dict = {}
    out_features = []
    for classy_gene in bog:
        gff_dict, chrom_len_dict = create_or_append_gff_feature(classy_gene,gff_dict,chrom_len_dict)
    for chrom in gff_dict:
        rec = gff_dict[chrom][0]
        print(len(rec))
        feature_list = gff_dict[chrom][1:]
        rec.features = feature_list
        out_features.append(rec)
    with open (out_file, 'w') as out_handle:
        GFF.write(out_features,out_handle)


















if __name__ == '__main__':
    infile = '/home/ndh0004/code/coge_tools/data_out/51576_52024_qac_ks_sep20_bag_of_genes.tsv'
    bog = read_bag_of_genes_classy(infile)
    #chrom_linkages_through_parent(bog,'/home/ndh0004/code/coge_tools/data_out/'
    #                                  '51576_52024_qac_ks_sep20_bag_of_genes_filter.gexf')
   # chrom_linkages_by_list(bog,'/home/ndh0004/code/coge_tools/data_out/'
    #                                  '51576_52024_qac_ks_sep20_bag_of_genes_filter_AB.gexf',
    """                          [
                            "scaffold306_size1284463_obj",
                            "scaffold591_size730508_subseq_1:450970_obj",
                            "scaffold93_size1230182_subseq_1:987274_obj",
                            "Super-Scaffold_10",
                            "Super-Scaffold_105",
                            "Super-Scaffold_1173",
                            "Super-Scaffold_12",
                            "Super-Scaffold_128",
                            "Super-Scaffold_13",
                            "Super-Scaffold_146",
                            "Super-Scaffold_147",
                            "Super-Scaffold_158",
                            "Super-Scaffold_202",
                            "Super-Scaffold_210",
                            "Super-Scaffold_2412",
                            "Super-Scaffold_25",
                            "Super-Scaffold_253",
                            "Super-Scaffold_266",
                            "Super-Scaffold_2765",
                            "Super-Scaffold_291",
                            "Super-Scaffold_307",
                            "Super-Scaffold_311",
                            "Super-Scaffold_39",
                            "Super-Scaffold_50",
                            "Super-Scaffold_70",
                            "Super-Scaffold_76",
                            "Super-Scaffold_84",
                            "Super-Scaffold_95"

                            ]
)"""

    include = ["scaffold306_size1284463_obj",
                            #"scaffold591_size730508_subseq_1:450970_obj",
                            #"scaffold93_size1230182_subseq_1:987274_obj",
                            "Super-Scaffold_10",
                            #"Super-Scaffold_105",
                            "Super-Scaffold_1173",
                            "Super-Scaffold_12",
                            "Super-Scaffold_128",
                            #"Super-Scaffold_13",
                            #"Super-Scaffold_146",
                            #"Super-Scaffold_147",
                            #"Super-Scaffold_158",
                            #"Super-Scaffold_202",
                            #"Super-Scaffold_210",
                            #"Super-Scaffold_2412",
                            #"Super-Scaffold_25",
                            #"Super-Scaffold_253",
                            #"Super-Scaffold_266",
                            #"Super-Scaffold_2765",
                            #"Super-Scaffold_291",
                            #"Super-Scaffold_307",
                            #"Super-Scaffold_311",
                            #"Super-Scaffold_39",
                            #"Super-Scaffold_50",
                            #"Super-Scaffold_70",
                            #"Super-Scaffold_76",
                            #"Super-Scaffold_84",
                            "Super-Scaffold_95"]

    new_bag = []
    for gene in bog :
        chrom = gene.ctregion.split("\t")[0]
        if chrom in include:
            new_bag.append(gene)
    convert_gff(new_bag,'/home/ndh0004/code/coge_tools/data_out/51576_52024_qac_ks_sep20_bag_of_genes.gff')

