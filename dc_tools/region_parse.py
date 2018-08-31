import pybedtools
import DagChainerParse.py
def parse_self_self(handle):
    """

    :param handle: Takes a diag.chain that is self v self
                   this is from CoGE. We really expect that
                   that most of the quality control on your s
                   settings has been done by this point.
    :return: [dict, bed]
             dict: links all syntetic regions to their match
                   provides a unique ID for each set of matching
                   syntentic regions.
             bed: bed is start and stop of each syntetic region.
                  this will be used downstream to determine which
                  regions to compare to progenitor genome.
    """
    f = open(handle)
    genes_in_block = []
    percent_id = []
    data = f.read().split("\n")
    if len(data[-1]) == 0:  # could just add an rstrip("\n") instead
        data = data[:-1]
    for line in data:
        if line[0] == '#':
            if block != 'Null':
                block_num, score, org_a, org_b, orient, num_gene_pairs = line.split('\t')
            block_num = int(block_num.lstrip('#'))
            score = float(score)
            num_gene_pairs = int(num_gene_pairs)
            block = Syn(block_num, score, org_a, org_b, orient, num_gene_pairs)
            genes_in_block = []
            percent_id = []
        elif line[0] != '#' and block != 'Null':
            gene_hit = gene_parser(line)
            org_a_percent_id = float(gene_hit.org_a_region.split('||')[-1])
            percent_id.append(org_a_percent_id)
    f.close()


if __name__ == '__main__':
    handle = '/home/ndh0004/Dropbox/Code/CoGe_parse/test_dag.all.go.aligncoords'
    parse_self_self(handle)

