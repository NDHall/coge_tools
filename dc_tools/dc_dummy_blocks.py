"""
It is an absolute pain to create DAGChainer output files for testing by hand. So here is a basic script that does. It
is fairly basic, but saves some time.

"""


def body_produce_hit(pgenome,pchrom,pstart, pstop,pgene,ppid):
    bodyp = '{pg}_{pchrom}\t{pchrom}||{pstart}||{pstop}||{pgene}||-1||CDS||2047451327||29837||{ppid}'.format(
        pg=pgenome,
        pchrom=pchrom,
        pstart=pstart,
        pstop=pstop,
        pgene=pgene,
        ppid=ppid)
    return bodyp

def start_stop_link(start,stop):
    return '{start}\t{stop}'.format(
        start=start,
        stop=stop
    )

def full_hit(bodya,link,bodyp,bodyend='220217\t1.000000e-250\t50'):
    return '{a}\t{l}\t{p}\t{e}\n'.format(
        a=bodya,
        l=link,
        p=bodyp,
        e=bodyend
    )


def header_produce(number,a_genome,achrom,p_genome,pchrom):
    header = '#1\t1600.0\t{ag}_{achrom}\t{pg}_{pchrom}\tf\t{n}\n'.format(
        n=number,
        ag=a_genome,
        pg=p_genome,
        achrom=achrom,
        pchrom=pchrom

    )
    return header
def file_dc(infile):
    out =[]
    f = open(infile)
    doc = f.read()
    f.close()
    doc = doc.rstrip("\n").split("\n")
    for line in doc:
        line = line.split("\t")
        if line[0] == '##':
            pass
        elif line[0] == '#':
            print(line)
            agenome, achrom, pgenome, pchrom, number = line[1:]
            out.append(header_produce(number,agenome, achrom, pgenome,pchrom))
        elif len(line) >0 :
            print(line)
            agene, astart, astop,apid,pgene,pstart,pstop,ppid = line
            abody = body_produce_hit(agenome,achrom,astart,astop,agene,apid)
            pbody = body_produce_hit(pgenome, pchrom, pstart, pstop, pgene, ppid)
            link = start_stop_link(astart,astop)
            out.append(full_hit(abody,link,pbody))
    return out


out_list=file_dc('/home/ndh0004/code/coge_tools/dc_tools/dc_gen.txt')

out='/home/ndh0004/code/coge_tools/test_data/v1_dc.txt'
fout = open(out,'w')
fout.write("".join(out_list))
fout.close()




