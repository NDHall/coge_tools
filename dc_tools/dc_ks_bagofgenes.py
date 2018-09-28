"""
Here we are going to take parent vs child ks/ka values and use a pre-determined ks value as a cutoff.
Any genes with that cutoff will be considered fair game. Output will consist of link gene, target region call,
percent id, original gene_name for post-hoc testing.

determine a or b genome is ref.

Set cutoff. Cutoff will be used to make calls. But we will report all syntenic gene relationships.

Calls for unused will be reported as well in case they will be useful in other contexts such as for future...


parse all gene positions that match ref.
For each gene position produce file

call    bgene    aregion    bregion  percent-id mean_percent-id ks  meanks  agene   agene agene_region  bgene \
bgene-region num_above_cutoff, totalnum, cutoff,



"""


from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import Bio.Alphabet as Alphabet
from Bio.Alphabet import IUPAC


try:
    import dc_ks_ka_to_bed
except:
    import dc_tools.dc_ks_ka_to_bed as dc_ks_ka_to_bed
try:
    import dc_parse
except:
    import dc_tools.dc_parse as dc_parse
import numpy as np

#make numpy more stringent.
np.seterr(invalid='raise')

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


class Gene_in_bag():
    def __init__(self, call, pgene, total_cregion,total_pregion, perid,meanperid,
                 ks, meanks,ka,meanka,
                 cgene, cgene_region,
                 pgene_region,syn_len,
                 pgene_orient, cgene_orient):
        self.call = call
        self.pgene = pgene
        self.ctregion = total_cregion
        self.ptregion = total_pregion
        self.perid = perid
        self.meanperid = meanperid
        self.ks = ks
        self.meanks = meanks
        self.cgene = cgene
        self.cgene_region = cgene_region
        self.pgene_region = pgene_region
        self.ka=ka
        self.meanka=meanka
        self.perid_len=None
        self.ka_len=None
        self.ks_len=None
        self.syn_len=syn_len
        self.cgene_orient=cgene_orient
        self.pgene_orient=pgene_orient


def parse_gene_list(gene_list,pid_cutoff,ks_cutoff,syn_len_cutoff,bag_of_genes_dict, strict_ks,
                    call ,ref='b',cutoff=3.0):
    """

    :param gene_list: When the DAGChainer file is parsed it is split into a block(the header) and a gene_list(the body)
                      We can do this by reading the whole DAGChainer file into memory. then splitting based on `#` then
                      splitting resulting blocks by line breaks. in the case of ks files. gene list == body[1:] because
                      the body comes with its own header that needs discarded.
    :param pid_cutoff: The percent id cutoff for the mean pid of a syntenic block. I set the default in arparse object.
                       default == 90.0 for this program. This is because we are hunting a recent AB divergence.
    :param ks_cutoff:  maximum acceptable mean ks value for a syntenic block for a gene to be considered. This
                       is really a bit redundant after quota align. But it can be used as a filter. Though it mostly
                       succeeds in over-filtering the results, since ks and ka reports are not output for every gene.
                       Even if it is calculated for every gene.
    :param syn_len_cutoff: Minimum number of genes used a syntenic region for it to be used in calls.
    :param bag_of_genes_dict: It is  dict of dict {gene:{'pass':[Geneinbag,...j], 'fail':[Geneinbag,...j]}  genes that
                              pass all cutoff filters are added 'pass' else: gene is added in 'fail'.
    :param strict_ks:   Bool if true ks_cutoff is used to filter genes that are added bag_of_genes_dict{}.
                        default==False. I have found that calls get really sparse with ks
    :param call:        Default value of call variabe. Useful in case call_ab is true.
    :param ref:         This takes the parent variable. The defaut is b which corresponds the the b genome in
                        DAGChainer output or column 2 of the DAGChainer report. Some may prefer to run SynMap with
                        the AB ref in col 'a'. This can specified here.
    :param cutoff:     codeml and by extension CoGe report ks values 0-inf, even though they are based on the
                       Nei and Gojobori 1986 approach. Though this is not the exact model used proxmial to the output.
                       Adjustments are made for models of nucleotide substitions. I have limited the max ks to 3.0 which
                       is pretty high since without accounting for models of substition this would be max. We don't want
                       highly different sites for this comparison.
    :return:
    """

    baglet_of_genes = []
    syn_len = len(gene_list)
    if ref == 'a':
        parent_start = dc_parse.region_parser(gene_list[0].org_a_region)
        child_start = dc_parse.region_parser(gene_list[0].org_b_region)
        parent_stop = dc_parse.region_parser(gene_list[-1].org_a_region)
        child_stop = dc_parse.region_parser(gene_list[-1].org_b_region)

    else:
        child_start = dc_parse.region_parser(gene_list[0].org_a_region)
        parent_start = dc_parse.region_parser(gene_list[0].org_b_region)
        child_stop = dc_parse.region_parser(gene_list[-1].org_a_region)
        parent_stop = dc_parse.region_parser(gene_list[-1].org_b_region)
    if parent_start.start  > parent_stop.stop :
        pstart = parent_stop.stop -1
        pstop = parent_start.start
    else :
        pstart = parent_start.start -1
        pstop = parent_stop.stop
    if child_start.start > child_stop.stop :
        cstart = child_stop.stop -1
        cstop = child_start.start
    else :
        cstart = child_start.start -1
        cstop = child_stop.stop
    ks_vals = []
    ka_vals = []
    pid_vals = []
    for gene in gene_list:


        if ref == 'a' :
            parent_region = dc_parse.region_parser(gene.org_a_region)
            child_region = dc_parse.region_parser(gene.org_b_region)

        else:
            parent_region = dc_parse.region_parser(gene.org_b_region)
            child_region = dc_parse.region_parser(gene.org_a_region)


        pgene = parent_region.gene
        cgene_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=child_region.chrom,
                  start=child_region.start -1 ,
                  stop=child_region.stop)
        pgene_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=parent_region.chrom,
                  start=parent_region.start -1 ,
                  stop=parent_region.stop)
        perid = float(child_region.percent_id)
        if type(perid) is float:
            pid_vals.append(perid)
        ks = gene.ks
        ka = gene.ka
        #print(gene.ks, type(gene.ks))
        meanka = None
        if type(gene.ks) is float and \
            gene.ks <= cutoff:
            ks_vals.append(gene.ks)
        if type(gene.ka) is float and \
            gene.ka <= cutoff:
            ka_vals.append(gene.ka)
        cgene = child_region.gene
        ctot_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=child_region.chrom,
                  start=cstart ,
                  stop=cstop)
        ptot_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=parent_region.chrom,
                  start=pstart,
                  stop=pstop)
        classy_gene = Gene_in_bag(
            call=call,
            pgene=pgene,
            total_cregion=ctot_region,
            total_pregion=ptot_region,
            perid=perid,
            ks=ks,
            ka=ka,
            meanka=None,
            meanks=None,
            cgene=cgene,
            cgene_region=cgene_region,
            pgene_region=pgene_region,
            syn_len=syn_len,
            meanperid=0,
            cgene_orient=int(child_region.orientation),
            pgene_orient=int(parent_region.orientation))
        baglet_of_genes.append(classy_gene)

    ka_len = None
    perid_len = None
    #print('ks vals:',ks_vals)
    if len(ks_vals) > 0:
        meanks = np.mean(ks_vals)
        ks_len = len(ks_vals)
    else:
        meanks = None
        ks_len = None
    if len(ka_vals) > 0:
        meanka = np.mean(ka_vals)
        ka_len = len(ka_vals)
    if len(pid_vals) >0:
        meanperid = np.mean(pid_vals)
        perid_len = len(pid_vals)
    # update all values.
    # add to dictionary in dictionary classify as a high quality or low quality based on ks mean for block and
    # or mean percent id per gene or per block.
    for classy_gene in baglet_of_genes:
        #print(meanperid)
        classy_gene.meanka = meanka
        classy_gene.meanks = meanks
        classy_gene.meanperid = meanperid
        classy_gene.perid_len = perid_len
        classy_gene.ka_len = ka_len
        classy_gene.ks_len = ks_len
        # this is where we filter for gene comparisons.
        #print(classy_gene.meanperid)
        if classy_gene.pgene not in bag_of_genes_dict:
            bag_of_genes_dict[classy_gene.pgene] = {'pass':[], 'fail':[]}
        if strict_ks == True:
            if classy_gene.meanks is not None and \
               classy_gene.meanks <= ks_cutoff and \
               classy_gene.meanperid >= pid_cutoff and \
               classy_gene.syn_len >= syn_len_cutoff:
                    bag_of_genes_dict[classy_gene.pgene]['pass'].append(classy_gene)
            else:
                    bag_of_genes_dict[classy_gene.pgene]['fail'].append(classy_gene)
        else:
            # now we just call based on perid cutoff and syntenic block length.
            if classy_gene.meanperid >= pid_cutoff and \
               classy_gene.syn_len >= syn_len_cutoff:
                    bag_of_genes_dict[classy_gene.pgene]['pass'].append(classy_gene)
            else:
                    bag_of_genes_dict[classy_gene.pgene]['fail'].append(classy_gene)

    return bag_of_genes_dict




def create_or_append_gff_feature(classy_gene,gff_dict):
    """
    This function was orginally written to take  just a single GeneinBag object and gff_dict.
    :param classy_gene: This is class Gene in Bag object
    :param chrom_len_dict: dictionary for determing lenghth of sequence.
    :return:dictionaries that will be used for creating gff file.
    """
    ct_chrom, ct_start, ct_stop = classy_gene.ctregion.split('\t')

    seq = SpoofSeq("ACTG")
    rec = SeqRecord(seq, ct_chrom)
    qualifiers = {"source": "prediction", "score": 10.0, "other": [classy_gene.ptregion.replace('\t', '_')],
                  "ID": classy_gene.ctregion.replace('\t', '_')}

    # Since we are just looking at the region the strand here is set as 0. This could be modified to include the
    # orientation from block, if block gets parsed. However, for these qac.ks files, that we are using, that
    # information does not make it to this stage.

    top_feature = SeqFeature(FeatureLocation(int(ct_start) + 1, int(ct_stop), ),
                             strand=0,
                             type='region',
                             qualifiers=qualifiers)
    sub_qualifiers = {"source": "prediction", "call": classy_gene.call, "linking_gene": classy_gene.pgene}

    cg_chrom, cg_start, cg_stop = classy_gene.cgene_region.split("\t")
    top_feature.sub_features = [SeqFeature(FeatureLocation(int(cg_start) + 1, int(cg_stop)),
                                           type='gene',
                                           strand=classy_gene.cgene_orient,
                                           qualifiers=sub_qualifiers)]
    if ct_chrom not in gff_dict:
        gff_dict[ct_chrom] = [rec, top_feature]
    else:
        gff_dict[ct_chrom].append(top_feature)


    return [gff_dict]

def create_or_append_gff_feature_AB_linkage(classy_gene_n,classy_gene_u,gff_dict,chrom_len_dict):
    """
    This function was written to take guess work out of adding lengths to class Seq that gff derives chrom len from.
    It includes init of SpoofSeq("ACTG"), which is class constructed for the express purpose of overriding __len__
    in the Seq() class. This allows us to assign an accurate length without reading the whole sequence in, which may be
    large and unwieldy or not proximate to the analysis.

    :param classy_gene_u linked gene
    :param classy_gene_n: This is class Gene in Bag object
    :param gff_dict: dictionary for creating final gff once everything is processed.
    :param chrom_len_dict: dictionary for deterring length of sequence.
    :return:dictionaries that will be used for creating gff file.
    """
    classy_gene = classy_gene_n
    ct_chrom, ct_start, ct_stop = classy_gene.ctregion.split('\t')
    seq = SpoofSeq("ACTG")
    if chrom_len_dict is not None \
        and ct_chrom in chrom_len_dict:
        seq.length = chrom_len_dict[ct_chrom]
    rec = SeqRecord(seq, ct_chrom)
    qualifiers = {"Source": "prediction", "Score": 10.0, "Note":[ "linking region:{l}".format(
        l=classy_gene.ptregion.replace('\t', '_')),'linked_region:{l}'.format(
        l=classy_gene_u.ctregion.replace("\t",'_'))],
                  "ID": classy_gene.ctregion.replace('\t', '_')
                    }
    # Since we are just looking at the region the strand here is set as 0. This could be modified to include the
    # orientation from block, if block gets parsed. However, for these qac.ks files, that we are using, that
    # information does not make it to this stage.

    top_feature = SeqFeature(FeatureLocation(int(ct_start) + 1, int(ct_stop), ),
                             strand=0,
                             type='region',
                             qualifiers=qualifiers)
    sub_qualifiers = {"Source": "prediction", "Note" :["call:{c}".format(c=classy_gene.call) ,
                                                       "linking_gene:{l}".format(l= classy_gene.pgene),
                                                       "linking_gene_region:{l}".format(l=classy_gene.pgene_region),
                                                       "linked_gene{l}".format(l= classy_gene_u.pgene),
                                                        "linked_gene_region{l}".format(l=classy_gene_u.cgene_region)]}

    print(classy_gene.cgene_orient)
    cg_chrom, cg_start, cg_stop = classy_gene.cgene_region.split("\t")
    top_feature.sub_features = [SeqFeature(FeatureLocation(int(cg_start) + 1, int(cg_stop)),
                                           type='gene',
                                           strand=int(classy_gene.cgene_orient),
                                           qualifiers=sub_qualifiers)]
    if ct_chrom not in gff_dict:
        gff_dict[ct_chrom] = [rec, top_feature]
    else:
        gff_dict[ct_chrom].append(top_feature)


    return gff_dict



def convert_gff(bog,out_file):
    """
    This is a simpler version of the above create_or_append_gff_feature_AB_linkage(). It does not adjust length of the
    the output. It is a helpful template, but ended up requring too much guess work about sequence length.

    :param bog: Takes bagof genes. This could be reduced several pythonic ways. Or you could call gff features on the
                fly using create_or_append_gff_feature() to call features as the arise duing read through or qualifying.

    :param out_file: handle to write gff 2
    :return:
    """
    chrom_len_dict = {}
    gff_dict = {}
    out_features = []
    for classy_gene in bog:
        gff_dict, chrom_len_dict = create_or_append_gff_feature(classy_gene,gff_dict)
    for chrom in gff_dict:
        rec = gff_dict[chrom][0]
        print(len(rec))
        feature_list = gff_dict[chrom][1:]
        rec.features = feature_list
        out_features.append(rec)
    with open (out_file, 'w') as out_handle:
        GFF.write(out_features,out_handle)



def write_dict_to_out_file(gene_dict, outfile):
    """

    :param gene_dict: this is the dictionary the gene information has been stored in 2 bags of genes.
            {gene:{'pass':[Geneinbag,...j], 'fail':[Geneinbag,...j]} It is a dictionary of dictionaries.
            This function writes all genes regardless of passing
    :param outfile: The out file already open that the dictionary contents are written into.
             g.call,
                    g.pgene,
                    g.ctregion, # 3 columns since it is 3 tab delimited fields
                    g.ptregion, # 3 columns since it is 3 tab delimited fields
                    g.perid,
                    g.meanperid,
                    g.ks,
                    g.meanks,
                    g.cgene,
                    g.cgene_region, # 3 columns since it is 3 tab delimited fields
                    g.pgene_region, # 3 columns since it is 3 tab delimited fields
                    g.ka,
                    g.meanka,
                    g.perid_len,
                    g.ka_len,
                    g.ks_len,
                    g.syn_len,
                    g.cgene_orient, #int negative or postive
                    g.pgene_orient  #int negative or postive
    :return:
    """
    all_genes = gene_dict['pass']
    if len( gene_dict['fail']) >0 :
     [ all_genes.append(x) for x in gene_dict['fail']]
    for g in all_genes:
        g_attrs = [ str(x) for x in [
                    g.call,
                    g.pgene,
                    g.ctregion,
                    g.ptregion,
                    g.perid,
                    g.meanperid,
                    g.ks,
                    g.meanks,
                    g.cgene,
                    g.cgene_region,
                    g.pgene_region,
                    g.ka,
                    g.meanka,
                    g.perid_len,
                    g.ka_len,
                    g.ks_len,
                    g.syn_len,
                    g.cgene_orient,
                    g.pgene_orient
                    ]]
        outstring = '\t'.join(g_attrs)+'\n'
        outfile.write(outstring)




def get_chrom_lens(in_len):
    """
    This is so that we can produce gff objects using biopython without loading the whole sequence into memory. We
    can just reference the contig name and the previousy measured length. It reads an input file formatted as
    <length>\t<contig name>\n
    <length>\t<contig name>\n
    ...


    :param in_len: space delimited files 'Len(int) chrom(str)\n'
    :return: dictionary { chrom : Len}
    """
    ret = {}
    with open(in_len,'r') as handle:
        lines = handle.read().rstrip("\n").split("\n")
        for line in lines:
            length, chrom = line.split(" ")
            assert chrom not in ret,'it appears there is a duplicate name in this set of chroms and lens' \
                                    ' {name} has appeared before.'.format(name=chrom)
            ret[chrom] = int(length)
    return ret




def parse_ks(infile,pid_cutoff,ks_cutoff,syn_len_cutoff, out_file,
             call,parent='b',strict_ks=True, call_ab=True,qac=True,in_len=None):
    """

    :param infile: This the handle for opening DAGChainer output.
    :param pid_cutoff: Minimum percent identitty to accepth
    :param ks_cutoff:  Maximum ks to accept.
    :param syn_len_cutoff: Minimum number of genes required to include genes from a syntentic block.
    :param out_file: This is a prefix that will be used for writitng out files.
    :param call: String this is the default input if no call is made for a gene. It can be used to label
                 species when calls are not being made.
    :param parent: This is the col of the DAGChainer file to use as the reference 'a'==0  and 'b'==1
    :param strict_ks: Only report calls that pass ks cutoffs. This drastically reduces call numbers.
    :param call_ab: Bool, call ab if the run is being used to separate A and B syntentic blocks.
    :param qac:  Was quota align used in the process if so True. This will directly affect how file is parsed.
    :param in_len: Default == None. If not none it is the handle for a tab delimited file that has chrom lengths.
                   These are turned into a dict of values so that we can spoof length of sequence.
                   structure ==
                   <length>\t<chrom name>\n
                   <length>\t<chrom name>\n
                   ...

    :return: output 5 files
            1. a set of _ab_gene.tsv calls.
            2. gff of A genes as children of syntenic regions.
            3. gff of B genes as children of syntenic regions.
            4. gff of all genes that were eligble to be called.
            5. _bag_of_genes.tsv. contains all the information for each gene stored in gene in bag variable.
    """
    bag_of_genes_dict = {}
    chroms_lens = None
    f = open(infile,'r')
    fout_calls = open('{out_file}_abcalls.tsv'.format(
        out_file=out_file
    ),'w')
    fout_total = open('{out_file}_bag_of_genes.tsv'.format(
        out_file=out_file
    ),'w')
    ks_file = f.read()

    if type(in_len) == str:
        chrom_lens = get_chrom_lens(in_len)
        print("in_loop")

    if qac is False:
        ks_file = ks_file.rstrip("\n").split('#')[2::]
    else:
        ks_file = ks_file.rstrip("\n").replace('\n###','\n#')
        ks_file = ks_file.split("#")[2::]
    for block, body in zip( ks_file[0::2],ks_file[1::2]):
        # we are not using block for this, right now. May want it later.
        body = body.rstrip("\n").split("\n")[1:] # we are dropping the comment line with column names.
        gene_list = dc_ks_ka_to_bed.ks_body_parser(body)
        bag_of_genes_dict = parse_gene_list(gene_list,pid_cutoff,ks_cutoff,syn_len_cutoff, bag_of_genes_dict,strict_ks,
                                            call,
                                            ref=parent)
    #Now that we have parsed the whole file into a bag of genes we can start comparing genes that pass our cutoffs.

    A_gff_dict = {}
    B_gff_dict = {}
    gff_dict = {}
    for gene in bag_of_genes_dict :
        #print(len(bag_of_genes_dict[gene]['pass']),len(bag_of_genes_dict[gene]['fail']))

        if call_ab is True:
            if len(bag_of_genes_dict[gene]['pass']) == 2:
                n,u = bag_of_genes_dict[gene]['pass']
                if n.perid > u.perid:
                    n.call = 'A'
                    u.call = 'B'
                    A_gff_dict = create_or_append_gff_feature_AB_linkage(n, u, A_gff_dict, chrom_lens)
                    B_gff_dict = create_or_append_gff_feature_AB_linkage(u, n, B_gff_dict, chrom_lens)
                elif n.perid < u.perid :
                    n.call = 'B'
                    u.call = 'A'
                    B_gff_dict = create_or_append_gff_feature_AB_linkage(n, u, B_gff_dict, chrom_lens)
                    A_gff_dict = create_or_append_gff_feature_AB_linkage(u, n, A_gff_dict, chrom_lens)
                else:
                    n.call = 'ambig'
                    u.call = 'ambig'
                # write to call file
                gff_dict = create_or_append_gff_feature_AB_linkage(n,u,gff_dict,chrom_lens)
                gff_dict = create_or_append_gff_feature_AB_linkage(u, n, gff_dict, chrom_lens)
                for called_g in [n,u]:
                    if called_g.call != 'ambig':
                        region = called_g.ctregion.replace('\t','_')
                        outstring = '{r}\t{c}\t{g}\n'.format(
                        r=region,
                        c=called_g.call,
                        g=called_g.pgene
                    )
                        fout_calls.write(outstring)



            write_dict_to_out_file(bag_of_genes_dict[gene],fout_total)



        else:
            write_dict_to_out_file(bag_of_genes_dict[gene], fout_total)
    fout_total.close()
    fout_calls.close()

    for g_dict, label in zip([A_gff_dict,B_gff_dict,gff_dict], ['A','B','Total']):
        out_features = []
        out_bed = '{o}_{l}.gff3'.format( o=out_file,
                                         l=label)
        for chrom in g_dict:
            rec = g_dict[chrom][0]
            print(len(rec))
            feature_list = g_dict[chrom][1:]
            rec.features = feature_list
            out_features.append(rec)

        with open(out_bed, 'w') as out_handle:
            GFF.write(out_features, out_handle)








if __name__ == '__main__':
    infile = '/home/ndh0004/code/coge_tools/test_data/ks_ivc.ks'
    infile = '/home/ndh0004/Downloads/' \
             'ks_analysis/51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks'
    out_file = '/home/ndh0004/code/coge_tools/data_out/bog_sep28_qaccall'
    infile = '/home/ndh0004/code/coge_tools/data' \
             '/51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all' \
             '.go_D20_g10_A5.aligncoords.Dm0.ma1.qac2.1.50.gcoords.ks'
    pid_cutoff = 90.0
    ks_cutoff = 3.0
    syn_len_cutoff = 5
    strict_ks = False
    qac = True
    call="cor"
    in_len = '/home/ndh0004/code/coge_tools/data/Ecor_PR202_scaffolds.len'
    parse_ks(infile,pid_cutoff,ks_cutoff,syn_len_cutoff,out_file,in_len=in_len,call=call,strict_ks=strict_ks, qac=qac)



