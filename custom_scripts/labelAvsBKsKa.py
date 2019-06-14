"""
NDH
This is a one off script to read in yamls and create output files based on cd label

"""
import yaml
import pickle


def ccMlMain():
    # start with dictionary that give annontation including kegg pathway.
    f = open('/home/ndh0004/Documents/keggPthCor/gene_dictv2.pckl', 'rb')
    ab_dict = pickle.load(f)
    f.close()

    cc_yaml = '/home/ndh0004/Documents/corKsKa/CorVsCorTotalNGMLv2.yaml'
    cor_kska = '/home/ndh0004/Documents/corKsKa/CorVsCorTotalNGML.tsv'
    e = open('/home/ndh0004/Documents/corKsKa/CorVsCorTotalNGML.error', 'w')

    tk = open(cor_kska,'w')
    ccyl = open(cc_yaml, 'r')
    kska = yaml.load(ccyl)

    tk.write('MLds\tMLdn\tNGds\tNGdn'
             '\tcall0,\tscafType0\tSitaProt0\tSitaKegg0\KeggOrth0\tseq0'\
             '\tcall1,\tscafType1\tSitaProt1\tSitaKegg1\KeggOrth1'\
              '\tseq1\tid\n')
    for key in kska:
        call = []
        for seq in kska[key]['seqs']:
            print (kska[key]['seqs'])
            key_seq = seq.replace('|','_').replace('-mRNA-1','')
            if key_seq in ab_dict:
                if len(ab_dict[key_seq])== 2:
                    ab_dict[key_seq] += ['na','na','na']
                    call.append(ab_dict[key_seq])
                elif len(ab_dict[key_seq])== 3:
                    ab_dict[key_seq] += ['na','na']
                    call.append(ab_dict[key_seq])
                elif len(ab_dict[key_seq]) == 5:
                    call.append(ab_dict[key_seq])
                else:
                    assert False, 'malformed dict entry {k}:{d}'.format(k=key_seq,d=ab_dict[key_seq])
            else:
                call.append(['na','na','na','na','na'])

        if len(kska[key]['seqs']) == 2:
            outstring = '{mds}\t{mdn}\t{nds}\t{ndn}\t{call0}\t{seq0}\t{call1}\t{seq1}\t{k}\n'.format(
                mds=kska[key]['dS']['ML'],
                mdn=kska[key]['dN']['ML'],
                nds=kska[key]['dS']['NG'],
                ndn=kska[key]['dN']['NG'],
                call0='\t'.join(call[0]),
                call1='\t'.join(call[1]),
                seq0=kska[key]['seqs'][0],
                seq1=kska[key]['seqs'][1].replace('|','_').replace('-mRNA-1',''),
                k=key
            )
            tk.write(outstring)
        else:
            e.write('key: {k}\n'.format(k=key))
    tk.close()
    e.close()

def cothersMlMain( cc_yaml, cor_kska, e ):
    # start with dictionary that give annontation including kegg pathway.
    f = open('/home/ndh0004/Documents/keggPthCor/gene_dictv2.pckl', 'rb')
    ab_dict = pickle.load(f)
    f.close()

    tk = open(cor_kska,'w')
    ccyl = open(cc_yaml, 'r')
    kska = yaml.load(ccyl)

    e = open(e, 'w')
    tk.write('MLds\tMLdn\tNGds\tNGdn'
             '\tcall0\tscafType0\tSitaProt0\tSitaKegg0\tKeggOrth0\tseq0'\
             '\tcall1\tscafType1\tSitaProt1\tSitaKegg1t\KeggOrth1'\
              '\tseq1\tid\n')
    for key in kska:
        call = []
        for seq in kska[key]['seqs']:
            print (kska[key]['seqs'])
            key_seq = seq.replace('|','_').replace('-mRNA-1','')
            if key_seq in ab_dict:
                if len(ab_dict[key_seq])== 2:
                    ab_dict[key_seq] += ['na','na','na']
                    call.append(ab_dict[key_seq])
                elif len(ab_dict[key_seq])== 3:
                    ab_dict[key_seq] += ['na','na']
                    call.append(ab_dict[key_seq])
                elif len(ab_dict[key_seq]) == 5:
                    call.append(ab_dict[key_seq])
                else:
                    assert False, 'malformed dict entry {k}:{d}'.format(k=key_seq,d=ab_dict[key_seq])
            else:
                call.append(['na','na','na','na','na'])

        if len(kska[key]['seqs']) == 2:
            outstring = '{mds}\t{mdn}\t{nds}\t{ndn}\t{call0}\t{seq0}\t{call1}\t{seq1}\t{k}\n'.format(
                mds=kska[key]['dS']['ML'],
                mdn=kska[key]['dN']['ML'],
                nds=kska[key]['dS']['NG'],
                ndn=kska[key]['dN']['NG'],
                call0='\t'.join(call[0]),
                call1='\t'.join(call[1]),
                seq0=kska[key]['seqs'][0],
                seq1=kska[key]['seqs'][1].replace('|','_').replace('-mRNA-1',''),
                k=key
            )
            tk.write(outstring)
        else:
            e.write('key: {k}\n'.format(k=key))
    tk.close()
    e.close()

def runcother():
    #Oropetium Run
    cc_yaml = '/home/ndh0004/Documents/corKsKa/OroVsCorTotalNGMLv2.yaml'
    cor_kska = '/home/ndh0004/Documents/corKsKa/OroVsCorTotalNGMLv2.tsv'
    e = '/home/ndh0004/Documents/corKsKa/OroVsCorTotalNGMLv2.error'
    cothersMlMain(cc_yaml, cor_kska, e)
    #Indica Run
    cc_yaml = '/home/ndh0004/Documents/corKsKa/Ind30KVsCorTotalNGMLv2.yaml'
    cor_kska = cc_yaml.replace('.yaml','tsv')
    e = cc_yaml.replace('.yaml','error')
    cothersMlMain(cc_yaml, cor_kska, e)

if __name__ == '__main__' :
    #    ccMlMain()
    runcother()


