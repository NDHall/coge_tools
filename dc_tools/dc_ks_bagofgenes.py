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

import dc_ks_ka_to_bed
import dc_parse
import numpy as np

#make numpy more stringent.
np.seterr(invalid='raise')

class Gene_in_bag():
    def __init__(self, call, pgene, total_cregion,total_pregion, perid,meanperid,
                 ks, meanks,ka,meanka,
                 cgene, cgene_region,
                 pgene_region,syn_len):
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


def parse_gene_list(gene_list,pid_cutoff,ks_cutoff,syn_len_cutoff,bag_of_genes_dict, strict_ks,call ,ref='b',cutoff=3.0):
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
            meanperid=0
        )
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

def write_dict_to_out_file(gene_dict, outfile):
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
                    g.syn_len
                    ]]
        outstring = '\t'.join(g_attrs)+'\n'
        outfile.write(outstring)


def parse_ks(infile,pid_cutoff,ks_cutoff,syn_len_cutoff, out_file, call,parent='b',strict_ks=True, call_ab=True,qac=True):
    bag_of_genes_dict = {}
    f = open(infile,'r')
    fout_calls = open('{out_file}_abcalls.tsv'.format(
        out_file=out_file
    ),'w')
    fout_total = open('{out_file}_bag_of_genes.tsv'.format(
        out_file=out_file
    ),'w')
    ks_file = f.read()
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
    for gene in bag_of_genes_dict :
        #print(len(bag_of_genes_dict[gene]['pass']),len(bag_of_genes_dict[gene]['fail']))
        if call_ab is True:
            if len(bag_of_genes_dict[gene]['pass']) == 2:
                n,u = bag_of_genes_dict[gene]['pass']
                if n.perid > u.perid:
                    n.call = 'A'
                    u.call = 'B'
                elif n.perid < u.perid :
                    n.call = 'B'
                    u.call = 'A'
                else:
                    n.call = 'eq'
                    u.call = 'eq'
                # write to call file
                for called_g in [n,u]:
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







if __name__ == '__main__':
    infile = '/home/ndh0004/code/coge_tools/test_data/ks_ivc.ks'
    infile = '/home/ndh0004/Downloads/' \
             'ks_analysis/51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks'
    out_file = '/home/ndh0004/code/coge_tools/test_out/bog_sep18_qaccall'
    infile = '/home/ndh0004/code/coge_tools/data' \
             '/51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all' \
             '.go_D20_g10_A5.aligncoords.Dm0.ma1.qac2.1.50.gcoords.ks'
    pid_cutoff = 90.0
    ks_cutoff = 2.0
    syn_len_cutoff = 5
    strict_ks = False
    qac = True
    call="cor"
    parse_ks(infile,pid_cutoff,ks_cutoff,syn_len_cutoff,out_file,call=call,strict_ks=strict_ks, qac=qac)



