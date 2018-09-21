"""
DAGChainer output is fairly involved.
Here we are parsing the output. I have
chosen to parse the table so that
each chrom/contig becomes its own object.
Containing Syntetic regions. Each Syntentic region will
in turn contain the reported genes that support it.
This scheme is a sensible way to interate through the data.
Another scheme would be to convert it all to data frame, which
may be useful for some other applications.

"""
#ToDo: make qac file converter. qac file ends in 2 comment lines which confounds final parsing. Easier
# and safer to convert qac to parseable version rather than try to keep expanding parser.
import pybedtools
#import networkx as nx
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

class Chrom:

    def __init__(self, chrom):
        self.chrom = chrom
        self.syn = []
        self.num_syn = 1
        self.syn_blocks = []

class Syn:
    """
    class for Syntetic blocks
    """


    def __init__(self,block_num, score, org_a, org_b, orient, num_gene_pairs ):
        self.block_num = int(block_num)
        self.score = float(score)
        self.org_a = org_a
        self.org_b = org_b
        self.orient = orient
        self.num_gene_pairs = int(num_gene_pairs)
        self.genes = []
        self.avg_per_id = 'Null'

class Gene:

    def __init__(self,org_a_chrom, org_a_region, org_a_start, org_a_stop,
            org_b_chrom, org_b_region, org_b_start, org_b_stop,
            evalue, accumulative_score ):
        self.org_a_chrom = org_a_chrom
        self.org_a_region = org_a_region
        self.org_a_start = int(org_a_start)
        self.org_a_stop = int(org_a_stop)
        self.org_b_chrom = org_b_chrom
        self.org_b_region = org_b_region
        self.org_b_start = int(org_b_start)
        self.org_b_stop = int(org_b_stop)
        self.evalue = evalue
        self.accumulative_scorer = accumulative_score

class Region:
    def __init__(self, chrom, start,stop,gene,orientation,seq_type,score1,score2,percent_id):
        self.chrom = chrom
        self.start = int(start)
        self.stop = int(stop)
        self.gene = gene
        self.orientation = int(orientation)
        self.seq_type = seq_type
        self.score1 = int(score1)
        self.score2 = float(score2)
        self.percent_id = float(percent_id)

class Percent_id_region_link:
    def __init__(self,self_node_id,parent_node_id,percent_id,raw_values,overlap):
        self.self_node_id = self_node_id
        self.parent_node_id = [parent_node_id]
        self.percent_id = [float(percent_id)]
        self.raw_values = [raw_values]
        self.overlap = [overlap]
    # so that  we can check number of matches
    def __len__(self):
        return len(self.parent_node_id)



def gene_parser(line):
    line = line.split('\t')
    #print(line)
    org_a_chrom, org_a_region, org_a_start, org_a_stop, \
    org_b_chrom, org_b_region, org_b_start, org_b_stop, \
    evalue, accumulative_score = line
    gene_hit = Gene(org_a_chrom, org_a_region, org_a_start, org_a_stop,
                    org_b_chrom, org_b_region, org_b_start, org_b_stop,
                    evalue, accumulative_score)
    return gene_hit
def region_parser( region):
    """

    :param region: takes region out of diag chain file. Essentially
                   these are blast hits
    :return: returns class region.
    """
    #print(region)
    chrom, start, stop, gene, orientation, seq_type, score1, score2, percent_id = region.split('||')
    region = Region(chrom, start, stop, gene, orientation, seq_type, score1, score2, percent_id)
    assert region.start <region.stop, "blast hits out of order in region %s"%(region)
    return region
def syntenic_block_header_parse(line):
    """

    :param line: this should be a syntenic header line from the dc output with 7 regions.
    :return: Syn class object.
    """
    block_num, score, org_a, org_b, orient, num_gene_pairs = line.split('\t')
    block_num = int(block_num.lstrip('#'))
    score = float(score)
    num_gene_pairs = int(num_gene_pairs)
    block = Syn(block_num, score, org_a, org_b, orient, num_gene_pairs)
    return block


def create_pybed4(chrom, start,stop,id):
    bed_str = "{chrom}   {start}  {stop}  {id}_id".format(
        chrom=chrom,
        start=start,
        stop=stop,
        id=id)


    new_pybed = pybedtools.BedTool(bed_str, from_string=True)
    return new_pybed

def create_pybed3(chrom, start,stop):
    #first bed used for determining overlaps.
    bed_str = "{chrom}   {start}  {stop}".format(
        chrom=chrom,
        start=start,
        stop=stop)
    new_pybed = pybedtools.BedTool(bed_str, from_string=True)
    return new_pybed



def deprc_bed_checker(raw_bed_coordinates,pybed):
    """
    bed checker looks for overlap among syntenic regions
    :param raw_bed_coordinates: These are take from start of the first gene in a syn block and stop of the last gene
                                in a syn block, they have been checked with assert to make sure start is smaller than
                                local stop, and vice-versa.
    :param pybed: BedTools object. used for checking to make sure syn regions do  not overlap.

    :return:[new_pybed, bed_asid, overlap_tags,overlap]
        :new_pybed:     new bed object that contains coordinates from current syn block
        :bed_asid:      bed as id, the chrom names and coordinates can be combined to make a unique identifier. This is done
                        so that it can be referenced again and added to the dictionary of connections.
        :overlap_tags:  a list of containing bed_asid, of all overlapping regions.
        :overlap:       boolean True overlap detected, False, not detected.

    """
    # here we need to check to make sure orientation is correct.
    # we worked to make sure maximum range was take start coordinate is smallest
    #   stop coordinate is largest. this is done using an assert statment when
    # now we need make sure that the orientation is correct. start < stop
    if raw_bed_coordinates[1] < raw_bed_coordinates[2]:
        start = raw_bed_coordinates[1]
        stop = raw_bed_coordinates[2]
    else:
        start = raw_bed_coordinates[2]
        stop = raw_bed_coordinates[1]

    bed_asid = '{chrom}_{start}_{stop}'.format(
        chrom=raw_bed_coordinates[0],
        start=start,
        stop=stop)
    overlap_tags = []
    new_pybed = create_pybed3(raw_bed_coordinates[0],start,stop)
    if len(pybed)== 0:
        overlap = False
    else :
        test_pybed = create_pybed3(raw_bed_coordinates[0],start,stop)
        tested_pybed = test_pybed.intersect(pybed,
                                            wo=True,
                                            f=0.90)
        if len(tested_pybed) == 0: # no overlaps just add to bed file
            new_pybed = create_pybed3(raw_bed_coordinates[0],start,stop)
            overlap = False
        elif len(tested_pybed) >0 : #Syn blocks appearing multiple times need to have the same id for each appearance.
            overlap = True
            for bed in tested_pybed: # get bed_asid tags for each overlapping region. ( we could also filter for length)
                overlap_id = bed[3:6] #bed element can be called by column if in pure bed form.
                overlap_id ='_'.join(overlap_id)

                assert overlap_id != bed_asid,"trouble parsing overlap flag."
                if overlap_id not in overlap_tags:
                    overlap_tags.append(overlap_id)
        new_pybed = pybed.cat(test_pybed,postmerge=False)
        for bed in [tested_pybed, test_pybed]:
            try:
                bed.delete_temporary_history(ask=False)
            except:
                pass #

    if len(pybed) >0:
        try:
            pybed.delete_temporary_history(ask=False)
        except:
            pass

    return [new_pybed, bed_asid, overlap_tags,overlap]


def self_bed_maker(raw_bed_coordinates, pybed,added_coords):
    """

    :param raw_coordinates: list of coordinates for syntetic block
    :param pybed: BedTools object
    :return: modified BedTools object containing syntetic region.
    """

    if raw_bed_coordinates[1] < raw_bed_coordinates[2]:
        start = raw_bed_coordinates[1]-1 # I am working on the principle that these are blast coords being converted to
                                         # bed coords.
        stop = raw_bed_coordinates[2]
    else:
        start = raw_bed_coordinates[2]-1
        stop = raw_bed_coordinates[1]
    bed_asid = '{chrom}_{start}_{stop}'.format(
                chrom=raw_bed_coordinates[0],
                start=start,
                stop=stop
    )

    added_coords[0].append(bed_asid)
    if bed_asid not in added_coords[1]:
        added_coords[1].append(bed_asid) # to prevent from adding identical regions multiple times.
        if len(pybed) == 0: # for first record added.
            new_pybed = create_pybed3(raw_bed_coordinates[0],start,stop)
        else:
            pybed_to_add = create_pybed3(raw_bed_coordinates[0],start,stop)
            new_pybed = pybed.cat(pybed_to_add, postmerge=False)
        try:
            pybed.delete_temporary_history(ask=False)
        except:
            pass
    else:
        new_pybed = pybed # nothing was added move on.

    return [new_pybed, added_coords]

def self_bed_produce(raw_bed_line_a, raw_bed_line_b, region_a, region_b, pybed, added_coords):
    """

    :param raw_bed_line_a: parsed CoGe Line
    :param raw_bed_line_b: parsed CoGe Line
    :param region_a:
    :param region_b:
    :param pybed:
    :return:
    """

    raw_bed_line_a.append(region_a.stop)  # add the final coordinates for the syntetic block for a and b.
    raw_bed_line_b.append(region_b.stop)  # it is okay if these are out of order because we will check that


    """
    next let's insert bed checking:
    """
    pybed, added_coords = self_bed_maker(raw_bed_line_b, pybed, added_coords)
    pybed, added_coords = self_bed_maker(raw_bed_line_a, pybed, added_coords)
    #print("bed:", len(pybed), "counter*2", counter * 2)
    return [pybed,added_coords]

def parse_1st_chainer_out(stem):
    """

    :param stem: raw DAGchainer output
    :return: [pybed,counter,added_coords]
             pybed: a pyBedTools object
             counter: number of syn blocks added.
             added_coords: [[ordered list of synblocks],
                            [All uniq nodes]]
    """
    #stem = '/home/ndh0004/Dropbox/Code/CoGe_parse/test_dag.all.go.aligncoords'
    f = open(stem)
    block = 'Null'
    raw_bed_line_a = []
    raw_bed_line_b = []
    counter = 1
    gene_counter = 0
    pybed = ''
    added_coords =[[],[]]


    data = f.read().split("\n")
    if len(data[-1]) == 0:  # could just add an rstrip("\n") instead
        data = data[:-1]
    for line in data:
        if line[0] == '#':
            if block != 'Null':

                pybed,added_coords = self_bed_produce(raw_bed_line_a, raw_bed_line_b,
                                                      region_a, region_b, pybed, added_coords
                                                      )

                block_num, score, org_a, org_b, orient, num_gene_pairs = line.split('\t')
                block_num = int(block_num.lstrip('#'))
                score = float(score)
                num_gene_pairs = int(num_gene_pairs)

                block = Syn(block_num, score, org_a, org_b, orient, num_gene_pairs)
                gene_counter = 0  # counts gene per syn block
                counter += 1  # counts num syn blocks. Used for dict matching
            else:
                block_num, score, org_a, org_b, orient, num_gene_pairs = line.split('\t')
                block_num = int(block_num.lstrip('#'))
                score = float(score)
                num_gene_pairs = int(num_gene_pairs)
                block = Syn(block_num, score, org_a, org_b, orient, num_gene_pairs)

        elif line[0] != '#' and block != 'Null':
            gene_hit = gene_parser(line)
            region_a = str(gene_hit.org_a_region)
            region_b = str(gene_hit.org_b_region)
            region_a = region_parser(region_a)
            region_b = region_parser(region_b)

            if gene_counter == 0:
                raw_bed_line_a = [region_a.chrom, region_a.start]
                raw_bed_line_b = [region_b.chrom, region_b.start]

                gene_counter += 1


    assert line[0] != '#', 'Malformed file. This file type should not end in \'#\''
    pybed, added_coords = self_bed_produce(raw_bed_line_a, raw_bed_line_b,
                                           region_a, region_b, pybed,
                                           added_coords)
    f.close()
    return [pybed,counter,added_coords]


def create_graph(list_of_edges):
    syn_graph =nx.Graph()
    for a , b in zip(list_of_edges[0::2], list_of_edges[1::2]):
        if syn_graph.has_edge(a,b) is False:
            syn_graph.add_edge(a,b)
        #print(a,b)
    return syn_graph

def add_intersect_to_syn_graph(syn_graph,intersect_bed):
    for overlap in intersect_bed:
        a = '{chrom}_{start}_{stop}'.format(
            chrom=overlap[0],
            start=overlap[1],
            stop=overlap[2]
        )
        b = '{chrom}_{start}_{stop}'.format(
            chrom=overlap[3],
            start=overlap[4],
            stop=overlap[5]
        )

        if a != b and syn_graph.has_edge(a,b) is False:
            syn_graph.add_edge(a,b)
    return syn_graph

def make_linkages(pybed,num_matches,added_coords,percent_overlap=0.70):
    total_added, unique_added = added_coords
    print( len(total_added))
    print("{matches} edges reported from from CoGe\n {unique} unique nodes.".format(
        matches=num_matches,
        unique=len(unique_added)))
    #create base graph
    syn_graph = create_graph(total_added) # create graph CoGe connections.
    #now find overlapping regions within bed.
    intersect_bed = pybed.intersect(pybed,
                    f=percent_overlap,
                    r=True,
                    wo=True)
    syn_graph = add_intersect_to_syn_graph(syn_graph, intersect_bed)

    return syn_graph


def parent_vs_self_bed_maker(raw_bed_coordinates, pybed, added_coords, percent_id,genes):
    if raw_bed_coordinates[1] < raw_bed_coordinates[2]:
        start = raw_bed_coordinates[
                    1] - 1  # I am working on the principle that these are blast coords being converted to
        # bed coords.
        stop = raw_bed_coordinates[2]
    else:
        start = raw_bed_coordinates[2] - 1
        stop = raw_bed_coordinates[1]
    bed_asid = '{chrom}_{start}_{stop}'.format(
        chrom=raw_bed_coordinates[0],
        start=start,
        stop=stop)
    #create id which will be mean^i,j,k...n
    mean_percent_id = np.mean(percent_id)
    raw_numbers = ','.join(str(x) for x in percent_id)
    raw_genes = ','.join((x.replace(',','#') for x in genes ))
    mean_and_raw = '{mean}^{raw_list}^{raw_genes}'.format(
        mean=mean_percent_id,
        raw_list=raw_numbers,
        raw_genes=raw_genes)
    added_coords[0].append(bed_asid)
    if bed_asid not in added_coords[1]:
        added_coords[1].append(bed_asid)  # to prevent from adding identical regions multiple times.
        if len(pybed) == 0:  # for first record added.
            new_pybed = create_pybed4(raw_bed_coordinates[0], start, stop,mean_and_raw)

        else:
            pybed_to_add = create_pybed4(raw_bed_coordinates[0], start, stop, mean_and_raw)
            new_pybed = pybed.cat(pybed_to_add, postmerge=False)
        try:
            pybed.delete_temporary_history(ask=False)
        except:
            pass
    else:
        new_pybed = pybed  # nothing was added move on.

    return [new_pybed, added_coords]


def parent_vs_self_bed_produce(raw_bed_line_a, raw_bed_line_b,
                               region_a, region_b, pybed, added_coords,
                               percent_id,genes):

    #print('regiona',region_a,'region_b',region_b)
    raw_bed_line_a.append(region_a.stop)  # add the final coordinates for the syntetic block for a and b.
    raw_bed_line_b.append(region_b.stop)  # it is okay if these are out of order because we will check that


    """
    next let's insert bed checking:
    """
    pybed, added_coords = parent_vs_self_bed_maker(raw_bed_line_b, pybed, added_coords, percent_id, genes)
    pybed, added_coords = parent_vs_self_bed_maker(raw_bed_line_a, pybed, added_coords, percent_id, genes)
    #print("bed:", len(pybed), "counter*2", counter * 2)
    return [pybed,added_coords]


def parse_parent_vs_self_syn(stem, parent='b'):
    """

    :param stem: file with DAG_chainer_output
    :param parent: this tells the function which region a or b contains parent, this is important, because we want to do
                   a 1 to 1 gene comparison. Followed by chi square test on region.
    :return: returns [pybed, counter, coords.
    """
    # stem = '/home/ndh0004/Dropbox/Code/CoGe_parse/test_dag.all.go.aligncoords'
    f = open(stem)
    block = 'Null'
    raw_bed_line_a = []
    raw_bed_line_b = []
    counter = 1
    gene_counter = 0
    pybed = ''
    added_coords = [[], []]
    percent_id =[]
    genes = []
    region_a = ''
    region_b = ''
    prev_line_start = ''


    data = f.read().split("\n")
    if len(data[-1]) == 0:  # could just add an rstrip("\n") instead
        data = data[:-1]
    for line in data:
        if line[0] == '#':
            if block != 'Null' and prev_line_start !='#':
                pybed, added_coords = parent_vs_self_bed_produce(raw_bed_line_a, raw_bed_line_b,
                                                       region_a, region_b, pybed, added_coords,
                                                       percent_id,genes)
                if len(line.split('\t')) == 7 :
                    block = syntenic_block_header_parse(line)
                else:
                    block = 'blank' # this value has to be set to != 'Null'

                gene_counter = 0  # counts gene per syn block
                counter += 1  # counts num syn blocks. Used for dict matching
            else:
                if len(line.split('\t')) == 7 :
                    block = syntenic_block_header_parse(line)
                else:
                    block = 'blank' # this value has to be set to != 'Null'

        elif line[0] != '#' and block != 'Null':
            gene_hit = gene_parser(line)
            region_a = str(gene_hit.org_a_region)
            region_b = str(gene_hit.org_b_region)
            region_a = region_parser(region_a)
            region_b = region_parser(region_b)
            percent_id.append(region_a.percent_id) # will be identical between the 2 hits
            if parent =='b':
                genes.append(region_b.gene)
            else:
                genes.append(region_a.gene)

            if gene_counter == 0:
                raw_bed_line_a = [region_a.chrom, region_a.start]
                raw_bed_line_b = [region_b.chrom, region_b.start]
                percent_id = [region_a.percent_id]
                if parent == 'b':
                    genes = [region_b.gene]
                else:
                    genes = [region_a.gene]
                gene_counter += 1
        prev_line_start = line[0]

    assert line[0] != '#', 'Malformed file. This file type should not end in \'#\''
    pybed, added_coords = parent_vs_self_bed_produce(raw_bed_line_a, raw_bed_line_b,
                                                       region_a, region_b, pybed, added_coords,
                                                       percent_id,genes)
    f.close()
    return [pybed, counter, added_coords]



def create_percent_id_region_link(line):
    """

    :param line: takes a line from a merge between bed4 and bed3
                 tricky bits are percent_id which is parsed here to extract total value of raw numbers and arth mean.
    :return: Percent_id_region_link
    """
    self_chrom, self_start, self_stop, parent_chrom, parent_start, parent_stop, percent_id, overlap = line
    self_node_name = '{chrom}_{start}_{stop}'.format(
        chrom=self_chrom,
        start=self_start,
        stop=self_stop)
    parent_node_name = '{chrom}_{start}_{stop}'.format(
        chrom=parent_chrom,
        start=parent_start,
        stop=parent_stop)

    percent, raw = percent_id.split('^')
    percent =float(percent)
    raw_nums = [float(x) for x in raw.rstrip("_id").split(",")]
    overlap = int(overlap)
    return Percent_id_region_link(self_node_name, parent_node_name, percent, raw_nums, overlap)



def create_node_percent_id_dict(percent_id_pybed, min_raw_overlap=2e4):
    node_dict = {}
    for line in percent_id_pybed:
        node = create_percent_id_region_link(line)

        if node.self_node_id not in node_dict and node.overlap[0] >= min_raw_overlap:
            node_dict[node.self_node_id] = node
        elif node.self_node_id in node_dict and node.overlap[0] >= min_raw_overlap:
            rec = node_dict[node.self_node_id]
            rec.parent_node_id.append(node.parent_node_id[0])
            rec.percent_id.append(node.percent_id[0])
            rec.raw_values.append(node.raw_values[0])
            rec.overlap.append(node.overlap[0])
    return node_dict






if __name__ == '__main__':
    import SynCompare

#Todo: add sort to bedtools
#Todo: change ref to smallest syn block
    self_vs_self_pybed,num_matches,added_coords = parse_self_vs_self_Syn(
        #'/home/ndh0004/Dropbox/Code/CoGe_parse/data/strict_test_self.dag_chainer'
        '/home/ndh0004/Dropbox/Code/CoGe_parse/SSacf_10.test' # for testing
        #'/home/ndh0004/Dropbox/Code/CoGe_parse/data/self_cor_51576_51576.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5 (1).aligncoords'
    )


    syn_graph = make_linkages(self_vs_self_pybed,num_matches,added_coords,percent_overlap=0.6 )
    parent_vs_self_pybed, counter, added_coords = parse_parent_vs_self_syn(
       # '/home/ndh0004/Dropbox/Code/CoGe_parse/data/tmp_self_v_out.dag_chainer'
       '/home/ndh0004/Dropbox/Code/CoGe_parse/SScaf_indica_v_cor.test' #for testing
       # '/home/ndh0004/Dropbox/Code/CoGe_parse/data/cor_vs_indica_51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords'
       )# currently a high overlap for testing.
    self_vs_self_pybed = self_vs_self_pybed.sort()
    parent_vs_self_pybed = parent_vs_self_pybed.sort()
    percent_id_pybed = self_vs_self_pybed.intersect(parent_vs_self_pybed
                                             ,wo=True
                                             ,f=0.20
                                              ,r=True
                                                    )
    #print("{matches} matches between parent_vs_self and self_vs self".format(matches=len(percent_id_pybed)))
    #print(percent_id_pybed[0][0]) # shows that we can in fact parse these quickly.

    node_dict = create_node_percent_id_dict(percent_id_pybed)



  #  for node in node_dict:
  #     print (node, node_dict[node].parent_node_id)




# next work through network.
    print( nx.connected_components(syn_graph))
    counter = 0
    f = open("/home/ndh0004/Dropbox/Code/CoGe_parse/tmp_v5.bed",'w') # for testing
    f.write('')
    f.close()
    del f
    f = open("/home/ndh0004/Dropbox/Code/CoGe_parse/tmp_v5.bed",'a')

    nx.nx.draw(syn_graph, with_labels=True)
    for con in nx.connected_components(syn_graph):


        counter += 1
        print("block_{index}".format(index=counter), len(con))
        out_bed = []
        if len(con) == 2:
            out_bed = SynCompare.match_dict_to_node(con, node_dict)
            if len(out_bed) == 2: # 2 connections syntenic
                un, nu = out_bed
                if len(un) == 1 and len(nu)== 1:
                    #ToDo: this block need turned into a function.
                    stat, p_val = ttest_ind(un.raw_values[0], nu.raw_values[0],
                                        equal_var=False)
                    print(un.parent_node_id,nu.parent_node_id,un.percent_id, nu.percent_id,stat,p_val)
                    if nu.percent_id[0] < un.percent_id[0]:
                        a_value = un.percent_id[0]
                    else:
                        a_value = nu.percent_id[0]
                    if p_val <=0.001 :
                        for x in [nu,un]:
                            SynCompare.write_linked_bed(f,x,a_value,counter,p_val)
                elif (len(un) > 1 or len(nu) >1): #it should be impossible to have an empty node object.
                    u_overlap, u_diff, u_raw_values = SynCompare.consolidated_raw_values(un)
                    n_overlap, n_diff,  n_raw_values = SynCompare.consolidated_raw_values(nu)
                    if u_overlap is False and n_overlap is False  and u_diff is False and n_diff is False:
                        stat, p_val = ttest_ind(u_raw_values, n_raw_values,
                                                equal_var=False)
                        print(un.parent_node_id, nu.parent_node_id, un.percent_id, nu.percent_id, stat, p_val)
                        if nu.percent_id[0] < un.percent_id[0]:
                            a_value = un.percent_id[0]
                        else:
                            a_value = nu.percent_id[0]
                        if p_val <= 0.001:
                            for x in [nu, un]:
                                SynCompare.write_linked_bed(f, x, a_value,counter)


    f.close()
    plt.show()











