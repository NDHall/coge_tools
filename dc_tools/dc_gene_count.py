import matplotlib.pyplot as plt
import numpy as np

def gene_count(line, blocks_counted, genes_counted, block, stop, counter, gene_counter):
    """

    :param line: takes line and determines if it is gene or new syn block. counts genes assoc.
                 with syn blocks.
    :param blocks_counted: list of blocks counted
    :param genes_counted: parallel list of genes per block counted.
    :param block: this tells program if its parsing from begining of file. If so leave as
                  'Null' else: set block manually.
    :return:
    """
    counter += 1

    if line[0] == '#':
        if block == 'Null':
            block = line.lstrip('#')
        elif block != 'Null' and gene_counter > 0:
            blocks_counted.append(block)
            genes_counted.append(gene_counter)
            gene_counter = 0
        else:
            block += line.lstrip('#')  # in case there are lines both starting with ## in a r
    elif block != 'Null' and stop == counter:
        blocks_counted.append(block)
        genes_counted.append(gene_counter+1) # because we

    elif line[0] != '#' and block != 'Null':
        gene_counter += 1
    return [blocks_counted, genes_counted, block, counter,gene_counter]


def gene_per_syn_count(stem):
    counter = 0
    gene_counter = 0
    block = 'Null'
    f = open(stem)
    data = f.read()
    data = data.rstrip('\n') # drop empty lines from end of file.
    data = data.split('\n')
    stop = len(data)
    genes_counted = []
    blocks_counted = []
    for line in data:
      blocks_counted, \
      genes_counted, \
      block, \
      counter,\
      gene_counter =  gene_count(line, blocks_counted, genes_counted, block,stop,counter, gene_counter)
    return [blocks_counted,genes_counted]

    #catch last line



def make_histogram(gene_nums,bins,out,show,dpi):
    plt.hist(gene_nums,bins=bins)
    if show != False:
        plt.show()
    plt.savefig(fname=out,
                dpi=dpi
                )


def message(string):
    if string == 'help':
        ret_string = """
        Takes a DAG_Chainer(dc) file and returns a historgram of number of genes in each syntenic block.
        """
    return ret_string


if __name__ == '__main__':
    import argparse
    import os

    cwd = str(os.getcwd())

    parser = argparse.ArgumentParser(description=message('help'))
    parser.add_argument('-i','--infile',
                        help='dc file',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-d','--dir',
                        help='directory to save results to, default is current working directory',
                        required=False,
                        default=cwd,
                        type=str,
                        dest='out_dir')
    parser.add_argument('--dpi',
                        help='dpi for output graph',
                        required=False,
                        default=300,
                        type=int,
                        dest='dpi')
    parser.add_argument('--bins',
                        help='size of bins for histogram',
                        required=False,
                        default=20,
                        type=int,
                        dest='bins')
    parser.add_argument('-s','--show',
                        help='show graph when made',
                        required=False,
                        default=False,
                        dest='show',
                        action="store_true")
    parser.add_argument('-o','--outfile',
                        required=True,
                        type=str,
                        dest="out")

    argv = parser.parse_args()

    out = '{dir}/{out}'.format(dir=argv.out_dir,out=argv.out.rstrip("/"))
    blocks, genes = gene_per_syn_count(argv.in_file)
    make_histogram(gene_nums=genes,bins=argv.bins,out=out,dpi=argv.dpi,show=argv.show)
