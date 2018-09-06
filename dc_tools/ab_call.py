"""
Script to parse and call AvsB sections of valid ab_syn_finder.py run.
"""
import scipy.stats as stats

def ab_call(in_file):
    region_call_tally = {}
    f = open(in_file)
    bed_genes = f.read().rstrip("\n").split("\n")
    f.close()
    for bed_gene in bed_genes:
        id, call, gene = bed_gene.split("\t")
        if id not in region_call_tally:
            region_call_tally[id] = {'A':0, 'B':0}
            region_call_tally[id][call] += 1
        else:
            region_call_tally[id][call] += 1


    return region_call_tally




def call_ab_region(region_out,region_call_tally):
    f_out = open(region_out, 'w')
    for region in region_call_tally:
        reg = region_call_tally[region]
        statistic, pvalue = stats.chisquare([reg['A'], reg['B']])
        call = 'ambig'
        if pvalue <= min_p and \
                reg['A'] + reg['B'] >= mincut and \
                reg['A'] != reg['B']:
            if reg['A'] > reg['B']:
                call = 'A'
            else:
                call = 'B'
        reg['call'] = call
        reg['pval'] = pvalue
        reg['stat'] = statistic
        f_out.write('{call}\t{chrom}\t{start}\t{stop}\t{A}\t{B}\t{csp}\t{css}\n'.format(
            chrom="_".join(region.split("_")[:-2]),
            start=region.split("_")[-2],
            stop=region.split("_")[-1],
            A=reg['A'],
            B=reg['B'],
            csp=pvalue,
            css=statistic,
            call=call))
    f_out.close()
    return region_call_tally

def add_region_to_chrom(id,region,chrom):
    if region['call'] == 'A' or \
        region['call'] == 'B':
        chrom['A'] += region['A']
        chrom['B'] += region['B']
        chrom['strong_call'] += 1
    else:
        chrom['a'] += region['A']
        chrom['b'] += region['B']
        chrom['ambig_call'] += 1
    chrom['region'].append(id)
    chrom['pvalue'].append(region['pval'])
    chrom['stat'].append(region['stat'])
    chrom['start'].append(int(id.split("_")[-2]))
    chrom['stop'].append(int(id.split("_")[-1]))
    chrom['call'].append(region['call'])



def create_chrom_tally(region_call_tally):
    chrom_call_tally= {}
    for id in region_call_tally:
        chrom = "_".join(id.split('_')[:-2])
        if chrom not in chrom_call_tally:
            chrom_call_tally[chrom] = {'A':0 ,'B':0,
                                       'a':0, 'b':0,
                                       'region':[],
                                       'call': [],
                                       'pvalue':[],
                                        'stat': [],
                                       'start' : [],
                                       'stop' : [],
                                       'strong_call':0,
                                       'ambig_call':0}
            add_region_to_chrom(id,region_call_tally[id],chrom_call_tally[chrom])
        else:
            add_region_to_chrom(id,region_call_tally[id],chrom_call_tally[chrom])
    return chrom_call_tally


def calculate_and_write(region_call_tally, out_pre):
    region_out = '{pre}_region_calls.tsv'.format(pre=out_pre)
    chrom_out = '{pre}_chrom_calls_all.tsv'.format(pre=out_pre)
    chrom_ambig = '{pre}_chrom_calls_ambig.tsv'.format(pre=out_pre)
    chrom_conf = '{pre}_chrom_calls_good.tsv'.format(pre=out_pre)
    chrom_prov = '{pre}_chrom_calls_prov.tsv'.format(pre=out_pre)
#write out region calls and assign call confidence
    region_call_tally = call_ab_region(region_out,region_call_tally)
#merge region calls into a chrom level calls.
    chrom_call_tally = create_chrom_tally(region_call_tally)
#create outfiles
    f_ambig = open(chrom_ambig,'w')
    f_all = open(chrom_out,'w')
    f_conf = open(chrom_conf, 'w')
    f_prov = open(chrom_prov,'w')
    header = '#call\tcall_qual\tnum_good_reg\tnum_ambig_reg\tchrom\tA_gene_calls\tB_gene_calls\tambig_a\tambig_b\n'
    for out_file in [f_ambig,f_all,f_conf,f_prov]:
        out_file.write(header)
#iterate over chrom calls and write results
    #firs assign calls
    for chrom in chrom_call_tally :
        cr = chrom_call_tally[chrom]
        cr['designation'] = 'ambig' # we set this a default value until it is explicity called.
        cr['chrom_level_call'] = 'Null'
    #strong calls
        if 'B' not in cr['call'] and \
            'A' in cr['call'] and \
            cr['ambig_call'] == 0 :
            cr['designation'] = 'strong'
            cr['chrom_level_call'] = 'A'
        elif 'A' not in cr['call'] and \
             'B' in cr['call'] and \
            cr['ambig_call'] == 0 :
            cr['designation'] = 'strong'
            cr['chrom_level_call'] = 'B'
        elif 'A' in cr['call'] and \
             'B' in cr['call'] and \
            cr['ambig_call'] == 0 :
            cr['designation'] = 'strong'
            cr['chrom_level_call'] = 'AB'
    #provisional calls
        #these calls could be worked with to make them more or less strict. At this point I am taking a broad view
        # of provisional calls.
        elif 'B' not in cr['call'] and \
                'A' in cr['call'] and \
                cr['ambig_call'] != 0:
            cr['designation'] = 'provisional'
            cr['chrom_level_call'] = 'A'
        elif 'A' not in cr['call'] and \
                'B' in cr['call'] and \
                cr['ambig_call'] != 0:
            cr['designation'] = 'provisional'
            cr['chrom_level_call'] = 'B'
        elif 'A' in cr['call'] and \
                'B' in cr['call'] and \
                cr['ambig_call'] != 0:
            cr['designation'] = 'provisional'
            cr['chrom_level_call'] = 'AB'
    #ambig calls
        #in this case not one good call was made on the chrom which will in most cases be a super scaffold.
        else :
            cr['designation'] = 'ambig'
            cr['chrom_level_call'] = 'ambig'
    #create outstring for writing
        out_string = '{call}\t{des}\t{sc}\t{ac}\t{chrom}\t{A}\t{B}\t{a}\t{b}\n'.format(
            chrom=chrom,
            A=cr['A'],
            B=cr['B'],
            call=cr['chrom_level_call'],
            des=cr['designation'],
            a=cr['a'],
            b=cr['b'],
            ac=cr['ambig_call'],
            sc=cr['strong_call']
        )
    #the universal file
        f_all.write(out_string)
    #write outstring to file based on chrom call criteria.
        if  cr['designation'] == 'strong':
            f_conf.write(out_string)
        elif cr['designation'] == 'provisional':
            f_prov.write(out_string)
        elif cr['designation'] == 'ambig':
            f_ambig.write(out_string)
        else:
            assert cr['designation'] != 'Null' , 'a chrom has not been called\n{chrom}'.format(chrom=chrom)
    for out_file in [f_ambig,f_all,f_conf,f_prov]:
        out_file.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="""
    takes syn_block output from ab_synfinder.py and calls A vs B for each region. If regions fall outside of cutoffs for
    regions output ambig is designated. If this occurs at the chrom/contig level, ambiguous calls can be found in 
    _ambig.tsv""")
    parser.add_argument('-i', '--infile',
                        help='per gene call from ab_syn_finder.py',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-o', '--outfile',
                        help="""ab_call.py produces mulitple out files each with their own name. Here you provided the 
                        base of the name ab_call.py will use to write out files for this analysis.""",
                        required=True,
                        type=str,
                        dest="out_pre")
    parser.add_argument('-p', '--pval_cutoff',
                        help="""Minimum acceptable p-value for a per region call. In which chi square analysis is used
                        to determine A vs B call. default == 0.05""",
                        required=False,
                        default=0.05,
                        type=float,
                        dest="pmin")
    parser.add_argument('-c', '--cutoff',
                        help="""Minimum acceptable number of genes used to call AvsB. It is suggested that this value is
                        informed by examination of WGD signatures from parent genomes. default == 15""",
                        required=False,
                        default=15,
                        type=int,
                        dest="cutoff")

    argv = parser.parse_args()

    in_file = argv.in_file
    out = argv.out_pre
    mincut = argv.cutoff
    min_p = argv.pmin
    region = ab_call(in_file)
    calculate_and_write(region,out )
