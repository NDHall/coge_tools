"""
The purpose of this script is to
create a dictionary of mRNA objects that can be used to classify
edger output


1. convert uniprot ids to kegg calls
2.merge kegg calls and uniprot ids with AB calls.

"""
from bioservices import KEGG
import pickle

def getKeggConvDict():
    """
    one time use function to create a pickled dictionary we can query with names.
    :return:
    """

    s = KEGG()
    convDb = s.conv('sita','ncbi-proteinid')
    print(convDb)
    print(type(convDb))
    s.find(convDb,'XP_004951166.1')

def convertUniProtToKegg(uni2gene):
    """

    :param uni2gene: 2 column list of gene name and best uniprot hit.
    :return: dictionary of {geneName: kegg}
    """

    u = UniProt()
    gene_to_uprot = {}
    f = open(uni2gene,'r')
    geneList = f.read().rstrip('\n').split('\n')
    for line in geneList :
        gene, uniprot = line.split('\t')
        assert gene not in gene_to_uprot,"Duplicate gene found {g}".format(g=gene)
        gene_to_uprot[gene] = uniprot
        print(u.mapping(fr='ID', to='KEGG_ID', query=""))
    f.close()



def main(uni2gene):
    convertUniProtToKegg(uni2gene)




if __name__ == '__main__':
    uni2gene = "/home/ndh0004/Documents/keggPthCor/testLink.txt"
    # main(uni2gene)
    getKeggConvDict()