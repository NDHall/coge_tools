#ToDo: create new ks/ka syn class and region class
try:
    import dc_parse
except:
    import dc_tools.dc_parse as dc_parse
import numpy as np
class Ks_syn(dc_parse.Syn):

    def process_float(self,val):
        indicator = False
        try :
            val = float(val)
            if type(val) == float:
                indicator = True
        except:
            val = str(val)
        return [val, indicator]



    def __init__(self,block_num, score, org_a, org_b, orient, num_gene_pairs,ks,ka ):
        super().__init__(block_num, score, org_a, org_b, orient, num_gene_pairs)
        self.ks, self.ks_is_float = self.process_float(ks)
        self.ka, self.ka_is_float = self.process_float(ka)


class Ks_gene(dc_parse.Gene):

    def process_float(self,val):
        indicator = False
        try :
            val = float(val)
            if type(val) == float:
                indicator = True
        except:
            val = str(val)
        return [val, indicator]


    def __init__(self, ks, ka, org_a_chrom, org_a_region, org_a_start, org_a_stop,
            org_b_chrom, org_b_region, org_b_start, org_b_stop,
            evalue, accumulative_score ):
        super().__init__(org_a_chrom, org_a_region, org_a_start, org_a_stop,
                     org_b_chrom, org_b_region, org_b_start, org_b_stop,
                     evalue, accumulative_score)
        self.ks, self.ks_is_float = self.process_float(ks)
        self.ka, self.ka_is_float = self.process_float(ka)







def bin_values(raw_range, dict_of_values, val):
    """

    :param raw_range:list of values to bin by
    :param dict_of_values: this is the dictionary that will hold values once bin has been assigned.
    :param val: value to be binned
    :return: dict_of_values that contain binned results.
    """
    counter = 1
    cont = True
    if val>raw_range[-1]:
        dict_of_values = add_value(dict_of_values, len(raw_range), val)
    else:
        for start,stop in zip(raw_range[:-1], raw_range[1:]):

            if counter == 1 and val < start:
                dict_of_values = add_value(dict_of_values,counter-1,val)
                cont = False
            elif counter != len(raw_range):
                if start <= val and val < stop :
                    #print("added")
                    dict_of_values = add_value(dict_of_values,counter,val)
                    cont = False
            elif counter == len(raw_range)-1:
                dict_of_values = add_value(dict_of_values,counter,val)
                cont = False

            counter += 1
    return dict_of_values


def ks_block_parser(ks_block):
    ks_block = ks_block.split('\t')
    block_num,score, org_a,org_b, orient, end_block1, end_block2 = ks_block
    end_block2 = end_block2.rstrip("\n").split(" ")
    ka = end_block2[-1]
    end_block1 = end_block1.split(" ")
    ks = end_block1[-1]
    num_gene_pairs = end_block1[0]
    ks_syn_block = Ks_syn( block_num, score, org_a, org_b, orient, num_gene_pairs, ks, ka )
    return ks_syn_block


def ks_body_parser(body):
    """
    Here we are parsing the body and changing into a list of Ks_gene objects.
    I want to keep this function so that it is portable. So I am willing to
    run through the list a second time to parse out desired values.
    :param body:
    :return:
    """
    return_list = []  #
    for gene in body:
        ks, ka, org_a_chrom, org_a_region, org_a_start, org_a_stop,\
        org_b_chrom, org_b_region, org_b_start, org_b_stop,\
        evalue, accumulative_score = gene.split("\t")
        gene_obj = Ks_gene(ks, ka, org_a_chrom, org_a_region, org_a_start, org_a_stop,
                                   org_b_chrom, org_b_region, org_b_start, org_b_stop,
                          evalue, accumulative_score)
        #print(gene_obj.ka_is_float,gene_obj.ka,gene_obj.ks,gene_obj.org_a_stop)
        return_list.append(gene_obj)
    return return_list



def region_to_bed(gene_list,counter, cutoff=2.0):
    ks = []
    ka = []
    bed_regions = []
    for gene in gene_list:
        if  gene.ka_is_float and gene.ks_is_float \
            and gene.ks <= cutoff and gene.ka <= cutoff :
            ks.append(gene.ks)
            ka.append(gene.ka)
    #print(gene_list)
    as1 = dc_parse.region_parser(gene_list[0].org_a_region) # regions are unchanged in ks/ka files.
    as2 = dc_parse.region_parser(gene_list[-1].org_a_region)
    bs1 = dc_parse.region_parser(gene_list[0].org_b_region) # regions are unchanged in ks/ka files.
    bs2 = dc_parse.region_parser(gene_list[-1].org_b_region)
    #print(ks,ka)
    for s1,s2 in zip([as1,bs1],[as2,bs2]):
        if s1.start > s2.start:
            start = s2.start -1
            stop = s1.stop
        else:
            start = s1.start - 1
            stop = s2.stop
        if len(ka) >0 and len(ks)>0:
            mean_ka = np.mean(ka)
            mean_ks = np.mean(ks)
        else:
            mean_ka = 'NA'
            mean_ks = 'NA'
        out_bed = '{chrom}\t{start}\t{stop}\t{pair_id}_id\t{num_genes}\t{num_measured}\t{mean_ks}\t{mean_ka}\t{raw_values_ks}\t{raw_values_ka}\n'.format(
            chrom=s1.chrom,
            start=start,
            stop=stop,
            num_genes=len(gene_list),
            num_measured=len(ka),
            mean_ka=mean_ka,
            mean_ks=mean_ks,
            raw_values_ka="^".join([str(x) for x in ka]),
            raw_values_ks = "^".join([str(x) for x in ks]),
            pair_id=counter
        )
        bed_regions.append(out_bed)
        #print(out_bed)
    return bed_regions








def ks_bed_produce(in_file,out,cutoff):
    """
    ks files here have 3 sections per syntenic block that we can extract using by stepping through the parsed
    file by 3 based on splitting on the '#'. There is one idiosyncrasy, the first line of file does not fit this pattern.
    By reading the whole file into memory and then splitting it will save some awkwardness. with final entries. These are
    typically large, files, but not so large that this approach will be untenable on a good laptop, or better yet
    a fatnode on local server.
    :return:
    """

    f = open(in_file,'r')
    ks_file = f.read()
    ks_file = ks_file.rstrip("\n").split('#')[2::] # we don't really want first line which is non-repeated comment.
    # step through by twos
    regions = [['#chom\tstart\tstop\tpair_id\tnum_genes\tnum_kska_genes\tmean_ks\tmean_ka\traw_ks\traw_ka\n'],]
    counter = 1
    for block, body in zip( ks_file[0::2],ks_file[1::2]):
        ks_syn_block = ks_block_parser(block)
        if ks_syn_block.ka_is_float == True : # right now concentrating on only blocks containing ks/ka values.
            body = body.rstrip("\n").split("\n")[1:] # we are dropping the comment line with column names.
            gene_list = ks_body_parser(body)
            #print(gene_list,body,block)
            regions.append(region_to_bed(gene_list,counter,cutoff=cutoff))


            #print(block)
        counter += 1
    fout = open(out,'w')
    for region in regions:
        fout.write("".join(region))
    fout.close()



#ToDo: make ks/ka bed producer.
#ToDo: ks/ka filter


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="""
    This script takes ks/ka file from CoGe tools and parsers it into bed regions. There is a cutoff parameter. To limit\
     the number ks/ka values. Since these can become unreasonably high. The mean ks/ka values are output in tsv form \
      for downstream parsing""")
    parser.add_argument('-i', '--infile',
                        help='dc file for ks ka comparisons',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-o', '--outfile',
                        required=True,
                        type=str,
                        dest="out")
    parser.add_argument('-c','--cut_off',
                        help='highest acceptable value for ks or ka. default == 2.0',
                        required=False,
                        default=2.0,
                        type=float,
                        dest='cut_off')
    argv = parser.parse_args()

    in_file = argv.in_file
    out = argv.out
    cutoff = argv.cut_off
    ks_bed_produce(in_file,out, cutoff)