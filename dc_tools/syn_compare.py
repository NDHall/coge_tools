import pybedtools
import DagChainerParse
from scipy.stats import ttest_ind
import numpy as np

def match_dict_to_node(con, node_dict):
    """

    :param con: subgroup of connections nodes are coordinates from self vs self
    :param node_dict: dictionary of information about each node obtained from self v out DAGchainer output
    :return: returns node information from dict if nodes are in dict and in con.
    """
    out_bed =[]
    for id in con:
        print(id)

        if id in node_dict:
            if True:
                out_bed.append(node_dict[id])

            else:
                pass
        else:
            pass
            print("{id} not found".format(id=id))
    return out_bed

def write_linked_bed(f,x,a_value,counter,p_val):
    if x.percent_id[0] == a_value:
        genome = 'A'
    else:
        genome = 'B'
    genome_out = {'{genome}^{counter}_id^{per_id}^{p_val}^{raw_values}'.format(
        genome=genome,
        counter=counter,
        per_id=x.percent_id[0],
        p_val=p_val,
        raw_values=x.raw_values[0]
    )}
    chrom = "_".join(x.self_node_id.split("_")[:-2])
    start = x.self_node_id.split("_")[-2]
    stop = x.self_node_id.split("_")[-1]
    f.write('{chrom}\t{start}\t{stop}\t{genome}\n'.format(
        chrom=chrom,
        start=start,
        stop=stop,
        genome=genome_out))

def consolidated_raw_values(node,expected_p_val=0.05,min_mean_sep=20.00):
    """

    :param node: takes a node containing multiple regions and determines if
                 they are overlapping using bedtools

    : param expected_p_val: min p value for test. All values greater expected_p_val will indicate that it is okay to
                            add lists together.
    :min_mean_sep:  minimum acceptable difference between list mean sizes. Additional protection against mismatch list
                    merging in case sample size is small.

    :return: [overlap: boolean, raw_values]
    """
    pybed = False
    overlap = False
    diff = False
    raw_values =[]
    for id, rw in zip(node.parent_node_id, node.raw_values):
        id = id.split("_")
        chrom = "_".join(id[:-2])
        start = str(id[-2])
        stop = str(id[-1])
        if len(raw_values) == 0:
            raw_values = raw_values + rw
        else :
            stat, p_val = DagChainerParse.ttest_ind(raw_values, rw,
                                    equal_var=False)

            if abs(np.mean(raw_values) - np.mean(rw)) >min_mean_sep \
                and p_val > expected_p_val:
                raw_values = raw_values + rw
            else :
                diff = True


        if pybed is False:
            pybed = DagChainerParse.create_pybed3(
                chrom, start,stop)
        else:
            new_pybed = DagChainerParse.create_pybed3(
                chrom, start,stop)
            intersect = len( pybed.intersect(new_pybed))
            if intersect >0 :
                overlap = True
    return [overlap, diff, raw_values]


