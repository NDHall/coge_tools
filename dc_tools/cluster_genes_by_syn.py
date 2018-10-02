"""
Here we are loading calls for AvsB back in from regions file and the intial bag of genes.
With these we will create a seed file. That can be grown by adding other DAGChainer files
to the process.

ho I put together the test data

(venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w Super-Scaffold_10 bog_sep28_qaccall_bag_of_genes.tsv | awk '{print $6}'| head  |  sort | uniq >tmp
(venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ less ../test_data_cluster/bog28_bag_of_genes.tsv
(venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f tmp bog_sep28_qaccall_bag_of_genes.tsv >../test_data_cluster/bog28_bag_of_genes.tsv




"""

import dc_tools.dc_ks_bagofgenes as bog
import viz.data_get as dg


class ReadRegion():
    def __init__(self, call, chrom, start,stop, a_gene_calls, b_gene_calls, chisq, chistat, passing):
        self.call = call
        self.chrom = chrom
        self.start=start
        self.stop = stop
        self.a_gene_calls = int(a_gene_calls)
        self.b_gene_calls = int(b_gene_calls)
        self.chisq = float(chisq)
        self.chistat = float(chistat)
        self.passing = passing
    def __str__(self):
        return '{c}\t{a}\t{o}'.format(
            c=self.chrom,
            a=self.start,
            o=self.stop

        )


def read_region(infile):
    """

    :param infile: This is the region file with no header produced by ab_call.py regions should only occur once per
                   per region. While their coordinates may overlap the genes used in calling the regions may not.
    :return: dictionary of objects class RegionRead()
    """
    f = open(infile,'r')
    read_f = f.read().rstrip("\n").split("\n")
    f.close()
    reg_dict = {}
    for line in read_f :
        call, \
        call_chrom, \
        start, \
        stop, \
        a_gene_calls, \
        b_gene_calls, \
        chisq, \
        chistat = line.split('\t')

        if call == 'A' or \
            call == 'B':
            passing = True
        else:
            passing = False

        rr = ReadRegion(call, call_chrom, start,stop, a_gene_calls, b_gene_calls, chisq, chistat, passing)
       # print(rr)
        assert rr not in reg_dict, """It appears that {r} has been repeated.""".format(r=rr)
        reg_dict[str(rr)]=rr
    return reg_dict





if __name__ == '__main__':
    infile = '/home/ndh0004/code/coge_tools/test_data_cluster/region_calls.tsv'
    bog_in = '/home/ndh0004/code/coge_tools/test_data_cluster/bog28_bag_of_genes.tsv'
    region_dict = read_region(infile)
    bog_dict = dg.read_bag_of_genes_classy(bog_in)
    # next read in bagof genes.
    # objectives
    #   produce dictionary of Eleusine regions
    #   {ER:{gene1:{region1:gene_a}}}
