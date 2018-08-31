def open_clean_for_append(out):
    """

    :param out: file to be opened for appending
    :return: returns clean version of file, so you don't accidentally keep appending to the same file.
    """

    f = open(out,'w')
    f.write('')
    f.close()
    del f
    f = open(out,'a')
    return f

def gene_filter(line, stop, counter, gene_counter,linesaver,cutoff,out,written):
    """

    :param line:
    :param stop:
    :param counter:
    :param gene_counter:
    :param linesaver:
    :param cutoff:
    :param out:
    :return:
    """
    counter += 1
    linesaver.append(line)
    if line[0] == '#':
        if len(linesaver)> 1 and linesaver[-2][0] == '#':
            pass # some Chainer files have 2 lines marking the start of each syn block.
        elif gene_counter >= cutoff:
            if written == 0 :
                out.write('\n'.join(linesaver[1:-1]))
            else:
                out.write('\n'.join(linesaver[:-1]) )
                print(linesaver)
                linesaver = ["\n"+linesaver[-1]]
                print(linesaver)
                gene_counter = 0
            written += 1
        elif counter == 1:
            pass
        else :
            linesaver = ["\n" + linesaver[-1]]
            gene_counter =0
    elif stop == counter and \
        gene_counter+1 > cutoff:
         out.write('\n'.join(linesaver) )
         out.write('\n')
    elif len(line)>0:
        gene_counter += 1

    return [linesaver,gene_counter, counter,written]


def gene_per_syn_filter(stem,out_stem,cutoff):
    linesaver = []
    counter = 0
    gene_counter = 0
    written = 0
    f = open(stem)
    out = open_clean_for_append(out_stem)
    data = f.read()
    data = data.rstrip('\n') # drop empty lines from end of file.
    data = data.split('\n')
    stop = len(data)
    for line in data:
        linesaver, gene_counter, counter, written =  gene_filter(line, stop, counter, gene_counter,linesaver,cutoff,out, written)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="""
    Reduce the syntenic blocks reported by raw number of genes contained in the syntenic block.
    Now CoGe synfinder reports genes in block that match. So this is just a post analysis way to 
    filter syntenic block size.
    """)
    parser.add_argument('-i','--infile',
                        help='dc file, tetraploid genome in a column and diploid genome in b column',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-o','--outfile',
                        help='outputfile name. No directory required.',
                        required=True,
                        type=str,
                        dest="out")
    parser.add_argument('-c','--cutoff',
                        help='minimum acceptable number of genes.',
                        required=True,
                        type=int,
                        dest="cutoff")
    argv = parser.parse_args()
    gene_per_syn_filter(argv.in_file, argv.out, cutoff=argv.cutoff)




