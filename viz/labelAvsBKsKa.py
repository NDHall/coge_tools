"""
NDH
This is a one off script to read in yamls and create output files based on cd label

"""
import yaml

def ngMain():
    kska_yaml = '/home/ndh0004/Documents/corKsKa/total_NG.yaml'
    #these go together for first run
    #a_list = '/home/ndh0004/Documents/corKsKa/A.list'
    #b_list = '/home/ndh0004/Documents/corKsKa/B.list'
    #total_kska = '/home/ndh0004/Documents/corKsKa/totalKsKa.tsv'
    a_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/A_all_calls.list'
    b_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/B_all_calls.list'
    total_kska = '/home/ndh0004/Documents/corKsKa/inferedCalls/ksKa_all_calls.tsv'

    #a_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/A_strong_calls.list'
    #b_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/B_strong_calls.list'
    #total_kska = '/home/ndh0004/Documents/corKsKa/inferedCalls/ksKa_strong_calls.tsv'



    tk = open(total_kska,'w')
    kyl = open(kska_yaml, 'r')
    al = open(a_list, 'r')
    bl = open(b_list, 'r')


    a = al.read().rstrip('\n').split('\n')
    b = bl.read().rstrip('\n').split('\n')

    kska = yaml.load(kyl)
    tk.write('ds\tdn\tcall\tseq0\tseq1\tid\n')
    for key in kska:
        call = []
        for seq in kska[key]['seqs']:
            if seq in a:
                call.append('a')
            elif seq in b:
                call.append('b')
            else:
                call.append('n')
        outstring = '{ds}\t{dn}\t{call}\t{seq0}\t{seq1}\t{k}\n'.format(
            ds=kska[key]['dS'],
            dn=kska[key]['dN'],
            call=''.join(call),
            seq0=kska[key]['seqs'][0],
            seq1=kska[key]['seqs'][1],
            k=key
        )
        tk.write(outstring)
    tk.close()
    kyl.close()
    bl.close()
    al.close()


def mlMain():
    kska_yaml = '/home/ndh0004/Documents/corKsKa/total_mlng.yaml'
    #these go together for first run
    #a_list = '/home/ndh0004/Documents/corKsKa/A.list'
    #b_list = '/home/ndh0004/Documents/corKsKa/B.list'
    #total_kska = '/home/ndh0004/Documents/corKsKa/totalKsKa.tsv'
    a_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/A_all_calls.list'
    b_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/B_all_calls.list'
    total_kska = '/home/ndh0004/Documents/corKsKa/inferedCalls/ksKa_all_calls.tsv'

    #a_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/A_strong_calls.list'
    #b_list = '/home/ndh0004/Documents/corKsKa/inferedCalls/B_strong_calls.list'
    #total_kska = '/home/ndh0004/Documents/corKsKa/inferedCalls/ksKa_strong_calls.tsv'



    tk = open(total_kska,'w')
    kyl = open(kska_yaml, 'r')
    al = open(a_list, 'r')
    bl = open(b_list, 'r')


    a = al.read().rstrip('\n').split('\n')
    b = bl.read().rstrip('\n').split('\n')

    kska = yaml.load(kyl)
    tk.write('MLds\tMLdn\tNGds\tNGdn\tcall\tseq0\tseq1\tid\n')
    for key in kska:
        call = []
        for seq in kska[key]['seqs']:
            if seq in a:
                call.append('a')
            elif seq in b:
                call.append('b')
            else:
                call.append('n')
        if len(kska[key]['seqs']) == 2:
            outstring = '{mds}\t{mdn}\t{nds}\t{ndn}\t{call}\t{seq0}\t{seq1}\t{k}\n'.format(
                mds=kska[key]['dS']['ML'],
                mdn=kska[key]['dN']['ML'],
                nds=kska[key]['dS']['NG'],
                ndn=kska[key]['dN']['NG'],
                call=''.join(call),
                seq0=kska[key]['seqs'][0],
                seq1=kska[key]['seqs'][1],
                k=key
            )
            tk.write(outstring)


if __name__ == '__main__' :
    mlMain()


