"""
Custom script
NDH
Purpose: to get list of unexpressed genes and to determine if blocks of non-expression exist.

"""



def groupByCluster():

    zero_indiv = []
    zero_contig_blocks = {}
    import pandas as pd
    last_hit = 0
    last_contig = 'Not a Contig'
    counts = '/home/ndh0004/Documents/cor_exp/raw_total_counts/gene.counts'
    df = pd.read_csv(counts,sep=" ",low_memory=False)
    print(len(df))
    contig_block = []
    for index, row in  df.iterrows() :
       row_sum = sum(row[2:])
       if row_sum == 0:

           zero_indiv.append(row['gene'])
           if (last_hit + 1 == index and last_contig == row['contig'] ) or \
                   (last_contig == 'Not a Contig' and last_hit == 0 ):
               contig_block.append(row['gene'])
           elif last_hit + 1 != index or last_contig != row['contig'] or index + 1 == len(df):
               assert len(contig_block) >0 , 'trouble with your row dude {c},{g}'.format(c=row['contig'], g=row['gene'])
               if last_contig not in zero_contig_blocks:
                   zero_contig_blocks[last_contig] = { len(contig_block) : [contig_block]}
               elif last_contig in zero_contig_blocks:
                   if len(contig_block) in zero_contig_blocks[last_contig]:
                       zero_contig_blocks[last_contig][len(contig_block)].append(contig_block)
                   else :
                          zero_contig_blocks[last_contig][len(contig_block)] = [contig_block]
               contig_block = [row['gene']]
               
           last_contig = row['contig']
           last_hit = index

    assert len(contig_block) >0 , 'trouble with your row dude {c},{g}'.format(c=row['contig'], g=row['gene'])
    if last_contig not in zero_contig_blocks:
        zero_contig_blocks[last_contig] = { len(contig_block) : [contig_block]}
    elif last_contig in zero_contig_blocks:
        if len(contig_block) in zero_contig_blocks[last_contig]:
            zero_contig_blocks[last_contig][len(contig_block)].append(contig_block)
        else :
               zero_contig_blocks[last_contig][len(contig_block)] = [contig_block]

    print(zero_contig_blocks[last_contig][1])
    f = open('/home/ndh0004/Documents/cor_exp/raw_total_counts/clusters.tsv','w')
    f.write('contig\tblock_len\tnum\n')
    for contig in zero_contig_blocks :
        for key in zero_contig_blocks[contig]:
            f.write('{c}\t{k}\t{l}\n'.format(c=contig, k=key, l=len(zero_contig_blocks[contig][key])))
    f.close()
    fout =open( '/home/ndh0004/Documents/cor_exp/raw_total_counts/single_zeros.tsv','w')
    for X in zero_indiv :
        fout.write(X+'\n')
    fout.close()

















if __name__ == '__main__':

    # Function to count up clusters on contigs.
    #    groupByCluster()

    #