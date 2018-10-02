import dc_tools.dc_ks_bagofgenes as bog
import pandas as pd


def read_bag_of_genes_classy(infile):
    bag_of_genes = []
    f = open(infile, 'r')
    raw_bag = f.read().rstrip('\n')
    for gene in raw_bag.split('\n') :
        gene = gene.split('\t')
        call, \
        pgene, \
        ctregion_chrom, \
        ctregion_start, \
        ctregion_stop, \
        ptregion_chrom, \
        ptregion_start, \
        ptregion_stop, \
        perid, \
        meanperid, \
        ks, \
        meanks, \
        cgene, \
        cgregion_chrom, \
        cgregion_start, \
        cgregion_stop, \
        pgregion_chrom, \
        pgregion_start, \
        pgregion_stop, \
        ka, \
        meanka, \
        perid_len, \
        ka_len, \
        ks_len, \
        syn_len,\
        cgene_orient,\
        pgene_orient   = gene

        classy_gene = bog.Gene_in_bag(call=call,
                                      pgene=pgene,
                                      total_cregion='{c}\t{sr}\t{sp}'.format(
                                          c=ctregion_chrom,
                                          sr=ctregion_start,
                                          sp=ctregion_stop
                                      ),
                                      total_pregion='{c}\t{sr}\t{sp}'.format(
                                          c=ptregion_chrom,
                                          sr=ptregion_start,
                                          sp=ptregion_stop
                                      ),
                                      perid=perid,
                                      meanperid=meanperid,
                                      ks=ks,
                                      meanks=meanks,
                                      cgene=cgene,
                                      cgene_region='{c}\t{sr}\t{sp}'.format(
                                          c=cgregion_chrom,
                                          sr=cgregion_start,
                                          sp=cgregion_stop),
                                      pgene_region='{c}\t{sr}\t{sp}'.format(
                                          c=pgregion_chrom,
                                          sr=pgregion_start,
                                          sp=pgregion_stop),
                                      ka=ka,
                                      meanka=meanka,
                                      syn_len=syn_len,
                                      pgene_orient=int(cgene_orient),
                                      cgene_orient=int(pgene_orient))
        classy_gene.perid_len = perid_len
        classy_gene.ka_len = ka_len
        classy_gene.ks_len = ks_len
        bag_of_genes.append(classy_gene)
    return bag_of_genes


def read_bog_to_pd(infile):
    ret = pd.read_csv(infile, sep="\t", header=None)
    ret.columns =  [
                    "call",
                    "pgene",
                    "ctregion_chrom",
                    "ctregion_start",
                    "ctregion_stop",
                    "ptregion_chrom",
                    "ptregion_start",
                    "ptregion_stop",
                    "perid",
                    "meanperid",
                    "ks",
                    "meanks",
                    "cgene",
                    "cgregion_chrom",
                    "cgregion_start",
                    "cgregion_stop",
                    "pgregion_chrom",
                    "pgregion_start",
                    "pgregion_stop",
                    "ka",
                    "meanka",
                    "perid_len",
                    "ka_len",
                    "ks_len",
                    "syn_len"]
    return ret

def read_abcall_to_pd(infile):
    ret = pd.read_csv(infile, sep="\t", header=None)
    ret.columns = ["region","pid","call","gene"]
    return ret


if __name__ == '__main__':
    print("""
    \t\tRunning data_get.py
    """)
