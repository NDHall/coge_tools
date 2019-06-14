# The raw names that included a pipe '|' needed to be fixed because they were dropped from
# the analysis.

ecorToSita_h = '/home/ndh0004/Documents/keggPthCor/Ecor2Sita.list'
ecorToSita_o = open(ecorToSita_h)
ecorToSita = ecorToSita_o.read().rstrip('\n').split('\n')

import pickle
f = open('/home/ndh0004/Documents/keggPthCor/gene_dict.pckl','rb')
ab_dict = pickle.load(f)
f.close()

#First lets load stragglers.

a_good = '/home/ndh0004/Documents/ABGenesByRegion/uniq_noPipeA_chrom_good.list'
b_good = '/home/ndh0004/Documents/ABGenesByRegion/uniq_noPipeB_chrom_good.list'
a_xo = '/home/ndh0004/Documents/ABGenesByRegion/uniq_noPipebog_sep28_qaccalled_merged_aba.list'
b_xo = '/home/ndh0004/Documents/ABGenesByRegion/uniq_noPipebog_sep28_qaccalled_merged_abb.list'


a_good_open = open(a_good)
a_good_list = a_good_open.read().rstrip('\n').split('\n')
a_good_open.close()
a_xo_open = open(a_xo)
a_xo_list = a_xo_open.read().rstrip('\n').split('\n')
a_xo_open.close()

b_good_open = open(b_good)
b_good_list = b_good_open.read().rstrip('\n').split('\n')
b_good_open.close()
b_xo_open = open(b_xo)
b_xo_list = b_xo_open.read().rstrip('\n').split('\n')
b_xo_open.close()





# add genes to ab_dict.

conflict_list = []
for gene_list,call,conf, in  zip([a_good_list, b_good_list, a_xo_list, b_xo_list],
                             ['a','b','a','b'],
                             ['good', 'good','xo','xo']):
    for gene in gene_list:
        if gene in ab_dict:
            assert ab_dict[gene][0] != call and ab_dict[gene][1] != conf , 'calls don\'t match {g},{dg}'.format(
                g=gene, dg=ab_dict[gene]            )

        else:
            ab_dict[gene] = [call, conf]
# Now add information about ecor
found = 0
found_list = []
missing = []
for line in ecorToSita :
    ecor, sita = line.split(' ')
    if ecor in ab_dict :
        if len(ab_dict[ecor]) == 2:
            ab_dict[ecor].append(sita)
            found += 1
            found_list.append(ecor)
    else:
        missing.append(ecor)
# load kegg module
from bioservices import KEGG
s = KEGG()
convDb = s.conv('sita','ncbi-proteinid')
convDb['ncbi-proteinid:YP_008815800']

# annotate kegg module
annotated =[]
no_joy_for_sita = []
for gene in found_list:
        sita = ab_dict[gene][-1]
        sita_q = 'ncbi-proteinid:{g}'.format(g=sita[:-2])
        if sita_q in convDb:
            ab_dict[gene].append(convDb[sita_q])
            annotated.append(gene)
        else:
            no_joy_for_sita.append(sita)
counter = 0

# get kegg pathways

for gene in annotated:
    counter += 1
    if (len(ab_dict[gene])) == 4:
        call, conf, sita, kSita = ab_dict[gene]
        keggObj = s.get(kSita)
        keggParse = s.parse(keggObj)
        ko = []
        if 'ORTHOLOGY' in keggParse.keys():
            ko = list(keggParse['ORTHOLOGY'].keys())
        else:
            ko = ['None']
        assert len(ko) == 1,'{ko:{k}\ngene:{g}'.format(ko=ko,g=gene)
        ab_dict[gene].append(ko[0])
    if (counter % 1000 ) == 0:
        print('Finshed:{c}'.format(c=counter))


f = open('/home/ndh0004/Documents/keggPthCor/gene_dictv2.pckl','wb')
pickle.dump(ab_dict,f)
f.close()

#success!