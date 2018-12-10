#!/usr/bin/env python
# coding: utf-8

#    

# # Making Dictionary to Classify genes

# * Load lists and create dictionary to populate.

# In[1]:


a_ambig = '/home/ndh0004/Documents/keggPthCor/a_ambig_genes.list'
a_good = '/home/ndh0004/Documents/keggPthCor/a_good_genes.list'
a_xo = '/home/ndh0004/Documents/keggPthCor/a_xo_genes.list'
b_ambig = '/home/ndh0004/Documents/keggPthCor/b_ambig_genes.list'
b_good = '/home/ndh0004/Documents/keggPthCor/b_good_genes.list'
b_xo = '/home/ndh0004/Documents/keggPthCor/b_xo_genes.list'
a_ambig_open = open(a_ambig)
a_ambig_list = a_ambig_open.read().rstrip('\n').split('\n')
a_ambig_open.close()
a_good_open = open(a_good)
a_good_list = a_good_open.read().rstrip('\n').split('\n')
a_good_open.close()
a_xo_open = open(a_xo)
a_xo_list = a_xo_open.read().rstrip('\n').split('\n')
a_xo_open.close()
b_ambig_open = open(b_ambig)
b_ambig_list = b_ambig_open.read().rstrip('\n').split('\n')
b_ambig_open.close()
b_good_open = open(b_good)
b_good_list = b_good_open.read().rstrip('\n').split('\n')
b_good_open.close()
b_xo_open = open(b_xo)
b_xo_list = b_xo_open.read().rstrip('\n').split('\n')
b_xo_open.close()
ab_dict = {}


# * Load dictionary

# In[2]:


conflict_list = []
for gene_list,call,conf, in  zip([a_good_list, b_good_list, a_xo_list, b_xo_list,a_ambig_list,b_ambig_list],
                             ['a','b','a','b','a','b'],
                             ['good', 'good','xo','xo','prov','prov']):
    for gene in gene_list:
        if gene in ab_dict:
            if ab_dict[gene][0] != call and                 ab_dict[gene][1] != conf :
                ab_dict[gene] += [call, conf]
                conflict_list.append(gene)
            else:
                #gtf prints out gene per exon, this can lead to inflated list.
                pass
        else:
            ab_dict[gene] = [call, conf]


# In[3]:


print(len(conflict_list))


# ## Add *Setaria italica* annontation
# 

# In[4]:


ecorToSita_h = '/home/ndh0004/Documents/keggPthCor/Ecor2Sita.list'
ecorToSita_o = open(ecorToSita_h)
ecorToSita = ecorToSita_o.read().rstrip('\n').split('\n')


# In[5]:


missing = []
found = 0
found_list =[]
for line in ecorToSita :
    ecor, sita = line.split(' ')
    if ecor in ab_dict:
        ab_dict[ecor].append(sita)
        found += 1
        found_list.append(ecor)
    else:
        missing.append(ecor)
print(found)


# In[6]:


len(missing)


# In[7]:


len(ab_dict)


# In[8]:


for missed in missing[0:20]:
    print(missed)


# In[9]:


short = []
for key in ab_dict:
    if len(ab_dict[key]) == 2:
                short.append(key)
print(len(short))


# Not all genes had hits that passed filtering step. Notice we have a fair number of masked genes. 

# In[10]:


print(short[0:10])
print('...')
print(short[-10:])


# * You can see this below. The filtering step was worth it to throw out low quality hits
# 
# * These missed genes are expected. 
# 
# egrep -w  PR202_augustus_masked-Super-Scaffold_16-processed-gene-3.823 Ecor_PR202_gene_names.fa protx_out_sita 
# in the Fasta->Ecor_PR202_gene_names.fa:>PR202_augustus_masked-Super-Scaffold_16-processed-gene-3.823-mRNA-1 gene=PR202_augustus_masked-Super-Scaffold_16-processed-gene-3.823 CDS=1-348
# In the Blast Results -> protx_out_sita:PR202_augustus_masked-Super-Scaffold_16-processed-gene-3.823-mRNA-1	XP_012701369.1	33.333	48	32	0	189	46	212	259	1.8	28.1
# protx_out_sita:PR202_augustus_masked-Super-Scaffold_16-processed-gene-3.823-mRNA-1	XP_004958466.1	40.625	32	19	0	123	28	191	222	3.2	27.3
# protx_out_sita:PR202_augustus_masked-Super-Scaffold_16-processed-gene-3.823-mRNA-1	XP_012704256.1	52.174	23	11	0	276	344	16	38	4.7	26.9
# 

# ## KEGG Pathways
# * Now let's add KO ids and names to each of these accessions. 

# In[12]:


#get a conversion dict
from bioservices import KEGG
s = KEGG()
convDb = s.conv('sita','ncbi-proteinid')
convDb['ncbi-proteinid:YP_008815800']


# In[13]:


no_joy_for_sita = []
annotated = []
for gene in found_list[0:10]:
    print(len(ab_dict[gene]))
    print(ab_dict[gene])
    

            


# In[14]:


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
            


# In[15]:


no_joy_for_sita[0:10]


# Not all genes are annontated in DAVID. Again, this is not a problem since we are only looking for a highly conserved pathway.

# In[16]:


len(annotated)


# In[18]:


print(annotated[0:10])


# In[19]:


counter = 0
for gene in annotated[0:10]:
    counter += 1
    print(counter,gene,ab_dict[gene])


# In[20]:


print(len(annotated))


# In[28]:


counter = 0
for gene in annotated[0:10]:
    counter += 1
    call, conf, sita, kSita = ab_dict[gene]
    keggObj = s.get(kSita)
    keggParse = s.parse(keggObj)
    ko = []
    if 'ORTHOLOGY' in keggParse.keys():
        ko = keggParse['ORTHOLOGY'].keys()
    else:
        ko = ['None']
    print(ko)
    print(counter,gene,ab_dict[gene])


# In[31]:


counter = 0
for gene in annotated[0:10]:
    counter += 1
    call, conf, sita, kSita = ab_dict[gene]
    keggObj = s.get(kSita)
    keggParse = s.parse(keggObj)
    ko = []
    if 'ORTHOLOGY' in keggParse.keys():
        ko = list(keggParse['ORTHOLOGY'].keys())
    else:
        ko = ['None']
    assert len(ko) == 1
    print(ko[0])


# In[33]:



counter = 0 
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


# In[34]:


for gene in annotated[0:10]:
    print(ab_dict[gene])


# * Save dictionary 

# In[37]:


import pickle
f = open('/home/ndh0004/Documents/keggPthCor/gene_dict.pckl','wb')
pickle.dump(ab_dict,f)
f.close()

