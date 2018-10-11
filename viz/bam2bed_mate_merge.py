"""
Simple script to merge bedfile coordinates.


"""
from pybedtools import BedTool

def main(handle):
    out_handle = handle.replace('.bed','merged.bed')
    fout = open(out_handle ,'w')
    read_merge = {}
    reads = BedTool(handle)
    for read in reads.sort():
        key = '/'.join(read[3].split('/')[:-1])
        if key not in read_merge :
            read_merge[key] = [read]
        else :
            read_merge[key].append(read)

    print("""
    Finshed reading Reads.\n
    {l} read_matches found""".format(l=len(read_merge)))

    for key in read_merge :
        if len(read_merge[key]) == 2 :
            n,u = read_merge[key]
            if n[0] == u[0]:
                if n.stop < u.start and \
                        n.stop < u.stop:
                    p5 = n
                    p3 = u
                elif  n.start > u.stop and \
                        n.start > u.start:
                    p5 = u
                    p3 = n
                else :
                    p5 = None
                    p3 = None
                    """
                    right now we are excluding all overlapping reads. 
                    """
                if p5 is not None :
                    start = p5.start
                    stop = p3.stop
                    assert start < stop , 'Malformed bed start: {start} stop:{stop}\n{bed1} : {bed2}'.format(
                                                                                 start=start,
                                                                                 stop=stop,
                                                                                 bed1=str(p5),
                                                                                 bed2=str(p3))

                    outstring = '{c}\t{start}\t{stop}\t{r}\t{score}\t.\n'.format(
                        c=p5[0],
                        start=p5[1],
                        stop=p3[2],
                        r=key,
                        score=p5[4])
                    fout.write(outstring)
    fout.close()








if __name__ == '__main__':
    for x in ['/home/ndh0004/code/coge_tools/data_out/MP_looseMapLocalNu30K.bed']:
        main(x)