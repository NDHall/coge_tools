"""
We are going to take region output from ab_calls.py and merge the output into bed regions. This will give us fine grain
resolution. Ultimately, it should not matter for homogeneous contigs, but it will be useful for AB contigs.
"""
import pybedtools
import dc_parse

def sep_valid_beds(lines):
    beds = {'A':[],'B':[], 'ambig':[], 'Achrom':[], 'Bchrom':[], 'ambigchrom' : []}

    for line in lines:
        call, chrom,start,stop,A,B,pval,stat = line.rstrip("n").split("\t")
        chromkey = '{call}chrom'.format(
            call=call        )
        beds[chromkey].append(chrom)
        bed_string = '{c}\t{start}\t{stop}'.format(
            c=chrom,
            start=start,
            stop=stop )
        if bed_string not in beds[call]:
            beds[call].append(bed_string)
    assert len(beds['A']) == len(beds['B']) , """
    Malformed file. File must contain equal number of A and B calls since to call an A you must have a B.
    """
    for X in beds:
        print (X)
    return beds


def create_merged_bedobject(beds,distance):
    bed_objs = {'Achrom':beds['Achrom'],'Bchrom':beds['Bchrom']}
    for call in ['A','B']:

        bed_string = "\n".join(beds[call])
        bed_objs[call] = pybedtools.BedTool(bed_string, from_string=True)
        bed_objs[call] = pybedtools.BedTool.merge(bed_objs[call].sort(),
                          d=distance)


    return bed_objs



def check_for_overlap(raw_beds,beds,distance):
    prev_vals =[]
    intersect_len = len(beds['A'].intersect(beds['B'], wo=True))
    while intersect_len == 0 and distance <=20000000:
        print(intersect_len, distance)
        prev_vals = [beds,distance]
        distance += 10000
        beds = create_merged_bedobject(raw_beds, distance)
        intersect_len = len(beds['A'].intersect(beds['B'], wo=True))
    return prev_vals



def sort_beds(beds,out_stem):
    A = 'Null'
    B = 'Null'
    ABa = 'Null'
    ABb ='Null'
    abb_out = '{out_stem}_abb.bed'.format(out_stem=out_stem)
    aba_out =  '{out_stem}_aba.bed'.format(out_stem=out_stem)
    a_out = '{out_stem}_a.bed'.format(out_stem=out_stem)
    b_out = '{out_stem}_b.bed'.format(out_stem=out_stem)
    for region in beds['A']:
        new_pybed = dc_parse.create_pybed3(str(region[0]), str(region[1]), str(region[2]))
        if str(region[0]) in beds['Bchrom']:
            #print(type(region))
            #print(type(AB),AB,region)
            if ABa == 'Null':
                ABa = new_pybed
            else :
                ABa = ABa.cat(new_pybed,output=aba_out)
        else:
            if A == 'Null':
                A = new_pybed
            else:
                A = A.cat(new_pybed,output=a_out)


    for region in beds['B']:
        new_pybed = dc_parse.create_pybed3(str(region[0]), str(region[1]), str(region[2]))
        if str(region[0]) in beds['Achrom']:
            if ABb == 'Null':
                ABb = new_pybed
            else:
                ABb = ABb.cat(new_pybed,output=abb_out)
        else:
            if B == 'Null' :
                B = new_pybed
            else:
                B = B.cat(new_pybed,output=b_out)
    print(len(ABb),len(ABa), len(B),len(A))
    aset = set(beds['Achrom'])
    bset = set(beds['Bchrom'])
    print(len(beds['A']),len(beds['B']))
    print(len(bset.intersection(aset)),len(aset),len(bset))










def parse_regions(in_file,out_stem):
    f = open(in_file)
    lines = f.read()
    lines = lines.rstrip("\n").split("\n")
    raw_beds = sep_valid_beds(lines)
    beds = create_merged_bedobject(raw_beds,600000)
    beds, distance = check_for_overlap(raw_beds, beds, 600000)
    sort_beds(beds,out_stem)





if __name__ == '__main__':
    parse_regions('/home/ndh0004/code/coge_tools/test_data/all_region_call.tsv',
                  '/home/ndh0004/code/coge_tools/test_out/bedtest_')
