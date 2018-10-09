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
        print(read.chrom)
        key = '/'.join(read[3].split('/')[:-1])
        if key not in read_merge :
            read_merge[key] = [read]
        else :
            read_merge[key].append(read)
    for key in read_merge :
        if len(read_merge[key]) == 2 :
            n,u = read_merge[key]
            print(n[0],u[0])
            if n[0] == u[0]:
                if n[1] < u[1] and \
                        n[2] <= u[1]:
                    p5 = n
                    p3 = u
                elif  n[1] > u[2] and \
                        n[1] > u[2]:
                    p5 = u
                    p3 = n
                else :
                    p5 = None
                    p3 = None
                    """
                    right now we are excluding all overlapping reads. 
                    """
                if p5 is not None :
                    outstring = '{c}\t{start}\t{stop}\t{r}\t{score}\t.\n'.format(
                        c=p5[0],
                        start=p5[1],
                        stop=p3[2],
                        r=key,
                        score=p5[4])
                    fout.write(outstring)
    fout.close()








if __name__ == '__main__':
    for x in ['/home/ndh0004/Downloads/tmp_passthrough_files/sort_by_name_MP_raw_DRR095902_Ecor_PR202_scaffoldsv2.bed']:
        main(x)