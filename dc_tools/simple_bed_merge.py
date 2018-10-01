"""
This is a simple script for merging bed objects by the 4th column.
"""

def main():
    infile= '/home/ndh0004/code/coge_tools/data_out/sorted_xo_chroms_oct1_2018.bed'
    genome_dict = '/home/ndh0004/code/coge_tools/data_out/cor_genome_sizes.txt'
    fg = open(genome_dict,'r')
    g = fg.read().rstrip('\n').split('\n')

    g_dict = {}
    for line in g :
        k,l = line.split('\t')
        g_dict[k]=l

    f = open(infile,'r')
    raw_bed = f.read().rstrip('\n').split('\n')

    # I am not dealing with headers here. This is a pretty specific case. This could easily be changed to exclude
    #headers, but I am not doing that right now.
    out_bed = []
    out_line = []
    current_call = None
    started_chroms = []
    for line in raw_bed:

        if line[0][0] == '#':
            #This provides the filtering for commented headers.
            pass
        else:
            chrom, start,stop,call = line.split('\t')
            if len(out_line) == 0:
                out_line = [chrom,'0']
                current_call = call
                last_stop = stop
                started_chroms.append(chrom)
            else :
                if current_call == call and chrom == out_line[0]:
                    last_stop = stop
                elif current_call != call and chrom == out_line[0]:
                    out_line += [last_stop,current_call]
                    out_bed.append(out_line)
                    current_call = call
                    out_line = [chrom,start]
                    last_stop = stop
                elif chrom != out_line[0]:
                    k = out_line[0]
                    l = g_dict[k]
                    out_line += [l,current_call] # want to call to the end of the chrom to simplify the bed in the next step.
                    out_bed.append(out_line)
                    current_call = call
                    out_line = [chrom,str(0)]
                    started_chroms.append(chrom)
                    last_stop=stop
    out_line += [last_stop, current_call]
    out_bed.append(out_line)

    outfile = '/home/ndh0004/code/coge_tools/data_out/pymerged_sorted_xo_chroms_oct1_2018.bed'
    fout = open(outfile,'w')
    for line in out_bed:
        fout.write("\t".join(line))
        fout.write("\n")
    f.close()
    fout.close()



if __name__ == "__main__":
    main()