""""
Looking for kegg pathways  in list of singleton genes.

"""

import pickle
# get list of kegg pathways
if __name__ == '__main__':
    f = open('/home/ndh0004/Documents/keggPthCor/gene_dictv2.pckl','rb')
    ab_dict = pickle.load(f)
    f.close()

    # load list of targeted singletons

    """ 
    had to adjust the names to exclude special characters that bash chokes on etc...
    ndh0004@IorekByrnison:~/Documents/keggPthCor/singletons$ egrep \> v4_labeled_interacting_singletons.fasta |\
     cut -b 2- | sed 's/[|:]/_/g' > v4_labeled_interacting_singletons.list

    """
    f= open('/home/ndh0004/Documents/keggPthCor/singletons/v4_labeled_interacting_singletons.list')
    siLi = f.read().rstrip('\n').split('\n')
    f.close()


    # load list of all singletons.
    for gene in siLi:
        gene = gene.replace('-mRNA-1','')
        if gene in ab_dict:
            print(gene, ' '.join(ab_dict[gene]))

        else:
            print(gene , ' '.join(['None'] * 5))

    # look at kegg pathways.




