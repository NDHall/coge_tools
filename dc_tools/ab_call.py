"""
Script to parse and call AvsB sections of valid ab_syn_finder.py run.
"""
import scipy.stats as stats

def ab_call(in_file):
    chrom_call_tally = {}
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

        chrom = "_".join(id.split('_')[:-2])
        if chrom not in chrom_call_tally :
            chrom_call_tally[chrom] = {'A' :0, 'B' : 0, 'region':[id]}
            chrom_call_tally[chrom][call] += 1
        else:
            chrom_call_tally[chrom][call] += 1
            chrom_call_tally[chrom]['region'].append(id)
    return [chrom_call_tally,region_call_tally]

def calculate_and_write(chrom_call_tally,region_call_tally, out_pre, mincut, min_p):
    region_out = '{pre}_region_calls.tsv'.format(pre=out_pre)
    chrom_out = '{pre}_chrom_calls_all.tsv'.format(pre=out_pre)
    chrom_ambig = '{pre}_chrom_calls_ambig.tsv'.format(pre=out_pre)
    chrom_conf = '{pre}_chrom_calls_conf.tsv'.format(pre=out_pre)
    f_out = open(region_out,'w')
    for region in region_call_tally:
        reg = region_call_tally[region]
        statistic, pvalue = stats.chisquare([reg['A'], reg['B']])
        call = 'ambig'
        if pvalue <= min_p and \
            reg['A']+reg['B'] >=mincut and \
            reg['A'] != reg['B']:
            if reg['A'] > reg['B']:
                call='A'
            else:
                call='B'
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
    f_ambig = open(chrom_ambig,'w')
    f_all = open(chrom_out,'w')
    f_conf = open(chrom_conf, 'w')

    for chrom in chrom_call_tally :
        cr = chrom_call_tally[chrom]
        fraction_a = cr['A']/(cr['A'] + cr['B'])
        fraction_b = cr['B']/(cr['A'] + cr['B'])

        equal = False
        if cr['A'] > cr['B'] :
            called = cr['A']
            uncalled = cr['B']
            called_str='A'
        elif cr['A'] < cr['B']:
            called = cr['B']
            uncalled = cr['A']
            called_str='B'
        else:
            equal = True
            called = 'Null'
            uncalled = 'Null'
            called_str='ambig'

        out_string = '{call}\t{chrom}\t{A}\t{B}\t{fa}\t{fb}\n'.format(
            chrom=chrom,
            A=cr['A'],
            B=cr['B'],
            fa=fraction_a,
            fb=fraction_b,
            call=called_str
        )
        f_all.write(out_string)

        if  uncalled <= mincut and \
            called >= mincut and \
            equal == False :

            f_conf.write(out_string)
        else:
            f_ambig.write(out_string)

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
    chrom, region = ab_call(in_file)
    calculate_and_write(chrom,region,out,mincut,min_p )