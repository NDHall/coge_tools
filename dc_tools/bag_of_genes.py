import dc_ks_bagofgenes

if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser(description="produces a bag of genes from DAGChainer output. This is useful "
                                                 "for creating  A vs B calls. By default 2 files are produced "
                                                 "<out_file>_abcalls.tsv and <out_file>_bag_of_genes.tsv "
                                                 "Both files preserve connections among"
                                                 "syntenic blocks through gene id of the genome designated parent. "
                                                 "By default, the parent genome is the second or 'b' genome entered "
                                                 "into SynMap CoGe. It can be changed to 'a' with flag --parent. ")
    parser.add_argument('--pid_cutoff',
                        help='minimum acceptable average for a syntenic block to be used in producing '
                             '<outfile>_ab_call.tsv. (default==90.0)',
                        required=False,
                        type=float,
                        default=90.0,
                        dest='pid_cutoff')
    parser.add_argument('--ks_cutoff',
                        help='minimum acceptable average ks average for a syntenic block to be used in producing '
                             'Note: that only values falling below 3.0 are averaged. This cutoff is set in the function '
                             ' parse_gene_list(gene_list ...cutoff=3.0). It can be modified directly in the function. '
                             '<outfile>_ab_call.tsv. (default==2.0)',
                        required=False,
                        default=2.0,
                        type=float,
                        dest='ks_cutoff')
    parser.add_argument('--syn_len_cutoff',
                        help='minimum number of genes in a syntenic block for it to be used in creating '
                             '<outfile>_ab_call.tsv. (default==5)',
                        required=False,
                        default=5,
                        type=int,
                        dest='syn_len_cutoff')
    parser.add_argument('--call',
                        help='sets the default value of the call column in the output. It will be useful when '
                             'looking at different syntentic blocks for which calls are not required.'
                             ' (default==\'N\')',
                        required=False,
                        default='N',
                        type=str,
                        dest='call')
    parser.add_argument('--parent',
                        help='this sets which column of the DAGChainer is to be used as the parent genome. DAGChainer '
                             'has an \'a\'(1) and a \'b\'(2) column.(default==\'b\') ' ,
                        required=False,
                        default='b',
                        type=str,
                        dest='parent')
    parser.add_argument('--qac_false',
                        help='use this flag if you produced DAGChainer output without using Quota Align in CoGe '
                             'There are slight differences in the delimiters.',
                        required=False,
                        default=True,
                        dest='qac',
                        action="store_true")
    parser.add_argument('--call_ab_false',
                        help='set flag if you don\'t want to produce <outfile>_ab_call.tsv',
                        required=False,
                        default=True,
                        dest='call_ab',
                        action="store_true")
    parser.add_argument('--strict_ks_true',
                        help='set strict_ks_true if you only want <outfile>_ab_call.tsv to contain ab calls passing '
                             'the ks filtering steps. Note: this really drops the number of calls.',
                        required=False,
                        default=False,
                        dest='strict_ks',
                        action="store_true")
    parser.add_argument('-o','--outfile',
                        required=True,
                        type=str,
                        dest="outfile")
    parser.add_argument('-i', '--infile',
                        help='ks/ka DAGChainer output from CoGe run with Quota Align, unless --qac_false. (see above)',
                        required=True,
                        type=str,
                        dest="infile")

    argv = parser.parse_args()

    dc_ks_bagofgenes.parse_ks(infile=argv.infile,
                              pid_cutoff=argv.pid_cutoff,
                              ks_cutoff=argv.ks_cutoff,
                              syn_len_cutoff=argv.syn_len_cutoff,
                              out_file=argv.outfile,
                              call=argv.call,
                              parent=argv.parent,
                              strict_ks=argv.strict_ks,
                              call_ab=argv.call_ab,
                              qac=argv.qac)