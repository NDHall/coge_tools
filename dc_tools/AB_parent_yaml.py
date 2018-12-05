import dc_tools.dc_ks_bagofgenes as bag
import viz.data_get as dg
import yaml

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


def create_genelink_dict(region_dict, bog):
    gene_dict = {'A':{}, 'B' :{}}
    for classy_gene in bog: # gene in bag == gib
        if classy_gene.call == 'A' or \
            classy_gene.call == 'B' :
            if classy_gene.ctregion in region_dict:  # not all regions available in bog are called for example singletons.
                classy_region = region_dict[classy_gene.ctregion]
                if classy_region.call not in gene_dict:
                    gene_dict[classy_region.call] = {}

                """
                Since we get syntenic blocks with fuzzy ends it is not possible to match parent syntenic region by exact coordinates
                alone. Here we do it by using the genes which are only called when unique to a gene block and occur twice as linking
                genes. This saves us using the coordinates and program like pybedtools. Which is excellent, but not appropriate for
                this particular problem. And not in line with the bag of gene approach we used to generate these outputs.
                """

                if classy_gene.cgene in gene_dict[classy_region.call]:
                    gene_dict[classy_region.call][classy_gene.cgene].append(classy_gene.pgene)
                    print('appendend classy gene')
                else:
                    gene_dict[classy_region.call][classy_gene.cgene] = [classy_gene.pgene]
                    print('added cg')

    return gene_dict


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
        print(rr)
        assert str(rr) not in reg_dict, """It appears that {r} has been repeated.""".format(r=rr)
        reg_dict[str(rr)]=rr
    return reg_dict


def main():
    infile = '/home/ndh0004/code/coge_tools/data_out/bog_sep28_qaccalled_region_calls.tsv'
    bog_in = '/home/ndh0004/code/coge_tools/data_out/bog_sep28_qaccall_bag_of_genes.tsv'
    region_dict = read_region(infile)
    bog = dg.read_bag_of_genes_classy(bog_in)
    gene_dict = create_genelink_dict(region_dict, bog)
    print(gene_dict)
    for key in gene_dict :
        print(gene_dict[key])
        out_handle = '{td}/{k}.yaml'.format(
            td='/home/ndh0004/Documents/gg_ortho_indica',
            k=key
        )
        yaml.dump(gene_dict[key],open(out_handle,'w'))








if __name__ == '__main__':
    main()