"""
Here we are loading calls for AvsB back in from regions file and the intial bag of genes.
With these we will create a seed file. That can be grown by adding other DAGChainer files
to the process.

ho I put together the test data

(venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w Super-Scaffold_10 bog_sep28_qaccall_bag_of_genes.tsv | awk '{print $6}'| head  |  sort | uniq >tmp
(venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ less ../test_data_cluster/bog28_bag_of_genes.tsv
(venv_btools) ndh0004@IorekByrnison:~/code/coge_tools/data_out$ egrep -w -f tmp bog_sep28_qaccall_bag_of_genes.tsv >../test_data_cluster/bog28_bag_of_genes.tsv




"""

import dc_tools.dc_ks_bagofgenes as bag
import viz.data_get as dg
import dc_tools.dc_ks_ka_to_bed as dkb



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


def parse_ks_bog_to_dict(infile, call, qac=True, pid_cutoff=50,ks_cutoff=3,syn_len_cutoff=5,strict_ks=False,parent='b'):
    bag_of_genes_dict = {}
    chroms_lens = None
    f = open(infile,'r')
    ks_file = f.read()

    if qac is False:
        ks_file = ks_file.rstrip("\n").split('#')[2::]
    else:
        ks_file = ks_file.rstrip("\n").replace('\n###','\n#')
        ks_file = ks_file.split("#")[2::]
    for block, body in zip( ks_file[0::2],ks_file[1::2]):
        # we are not using block for this, right now. May want it later.
        body = body.rstrip("\n").split("\n")[1:] # we are dropping the comment line with column names.
        gene_list = dkb.ks_body_parser(body)
        bag_of_genes_dict = bag.parse_gene_list(gene_list,
                                                pid_cutoff,
                                                ks_cutoff,
                                                syn_len_cutoff,
                                                bag_of_genes_dict,
                                                strict_ks,
                                                call,
                                                ref=parent,
                                                cutoff=3.0)
    return bag_of_genes_dict

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


def create_seed_dict(region_dict, bog):
    seed_dict = {}
    gene_dict = {}
    for classy_gene in bog: # gene in bag == gib
        if classy_gene.call == 'A' or \
            classy_gene.call == 'B' :
            if classy_gene.ctregion in region_dict:  # not all regions available in bog are called for example singletons.
                classy_region = region_dict[classy_gene.ctregion]

                """
                Since we get syntenic blocks with fuzzy ends it is not possible to match parent syntenic region by exact coordinates
                alone. Here we do it by using the genes which are only called when unique to a gene block and occur twice as linking
                genes. This saves us using the coordinates and program like pybedtools. Which is excellent, but not appropriate for
                this particular problem. And not in line with the bag of gene approach we used to generate these outputs.
                """

                if classy_gene.pgene in gene_dict:
                    ptregion = gene_dict[classy_gene.pgene]
                else:
                    ptregion = classy_gene.ptregion
                    gene_dict[classy_gene.pgene] = classy_gene.ptregion

                if ptregion not in seed_dict:
                    seed_dict[ptregion] = {classy_gene.pgene: {
                        str(classy_region): {'region': classy_region,
                                             'gene': classy_gene}
                    }}

                elif ptregion in seed_dict and \
                        classy_gene.pgene not in seed_dict[ptregion]:
                    seed_dict[ptregion][classy_gene.pgene] = {
                        str(classy_region): {'region': classy_region,
                                             'gene': classy_gene}
                    }
                elif ptregion in seed_dict and \
                        classy_gene.pgene in seed_dict[ptregion] and \
                        str(classy_region) not in seed_dict[ptregion][classy_gene.pgene]:
                    seed_dict[ptregion][classy_gene.pgene][str(classy_region)] = {
                        'region': classy_region,
                        'gene': classy_gene
                    }



                elif classy_gene.ptregion in seed_dict and \
                        classy_gene.pgene in seed_dict[classy_gene.ptregion] and \
                        str(classy_region) in seed_dict[classy_gene.ptregion][classy_gene.pgene]:
                    assert False, 'repeated gene for region that has been called. All genes should occur' \
                                  'only once. Repeated Region : {g}'.format(g=str(classy_region))
    return [seed_dict, gene_dict]

def add_dagchainer_output_to_seed_dict(seed_dict,gene_dict,add_dict):
    added_cgene = {}
    for key in add_dict :
        for classy_gene in add_dict[key]['pass']:
            if classy_gene.pgene in gene_dict and \
                classy_gene.cgene not in added_cgene:
                ptregion = gene_dict[classy_gene.pgene]
                # now we can add to seed_dict.
                seed_dict[ptregion][classy_gene.pgene][classy_gene.ctregion] = {
                    'region': 'uncalled',
                    'gene'  : classy_gene
                }
                added_cgene[classy_gene.cgene] = [[ptregion,classy_gene.pgene,classy_gene.ctregion,False]]
            elif classy_gene.pgene in gene_dict and \
                classy_gene.cgene in added_cgene:
                added_cgene[classy_gene.cgene].append([ptregion,classy_gene.pgene,classy_gene.ctregion,True])
                for entry in added_cgene[classy_gene.cgene]:
                    ptregion, pgene, ctregion,removed = entry
                    if removed == False :
                        seed_dict[ptregion][pgene][ctregion] = {
                    'region': 'uncalled',
                    'gene'  : 'removed'
                        }
                        entry[-1] = True
    return seed_dict


if __name__ == '__main__':
    infile = '/home/ndh0004/code/coge_tools/test_data_cluster/region_calls.tsv'
    bog_in = '/home/ndh0004/code/coge_tools/test_data_cluster/bog28_bag_of_genes.tsv'
    region_dict = read_region(infile)
    bog = dg.read_bag_of_genes_classy(bog_in)
    addfile = '/home/ndh0004/Downloads/tmp_passthrough_files/' \
              '28588_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac1.1.50.gcoords.ks'

    call = 'setv'
    """
    Now that we have bog and region_dict we can begin making matches. 
    We will begin by moving through the bog and building a new dict.
    """

    seed_dict, gene_dict = create_seed_dict(region_dict, bog)




    for call, addfile in zip(['setv','seti','oro'],['/home/ndh0004/Downloads/tmp_passthrough_files/28588_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac1.1.50.gcoords.ks',
'/home/ndh0004/Downloads/tmp_passthrough_files/28806_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac1.1.50.gcoords.ks',
'/home/ndh0004/Downloads/tmp_passthrough_files/'
'51527_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac3.1.50.gcoords.ks']) :

        add_dict = parse_ks_bog_to_dict(addfile, call)
        print('working on :{c}'.format(c=call), len(add_dict))

        seed_dict = add_dagchainer_output_to_seed_dict(seed_dict,gene_dict,add_dict)


    for pr in seed_dict : #pr= parent region:
        print(pr, len(seed_dict[pr]))
        for pg in seed_dict[pr]: #pg= parent gene
            out_list = [pg]
            for cg in seed_dict[pr][pg]:
                if seed_dict[pr][pg][cg]['region'] != 'uncalled':
                    call = seed_dict[pr][pg][cg]['region'].call
                else:
                    call = seed_dict[pr][pg][cg]['gene'].call
                out_list +=  ['{c}_{g}'.format( g=seed_dict[pr][pg][cg]['gene'].cgene,
                                                c=call)]
            print("\t".join(out_list))


    # objectives
    #   produce dictionary of Eleusine regions
    #   {ER:{gene1:{region1:gene_a}}}
