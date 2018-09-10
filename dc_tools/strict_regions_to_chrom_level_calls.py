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

def filter_bed(bed_obj, keep_bed, exclude_bed):
    new_bed = []
    bed_string = "\n".join(bed_obj)
    bed_obj = pybedtools.BedTool(bed_string, from_string=True)
    bed_obj = bed_obj.sort()
    if keep_bed is not None:
        #keep_bed = pybedtools.BedTool(keep_bed)
        bed_obj = bed_obj.intersect(keep_bed,
                                    wa=True)
    if exclude_bed is not None:
        bed_obj = bed_obj.intersect(exclude_bed,
                                    wa=True,
                                    v=True,
                                    r=0.90) # this is hard coded here, because it will be likely that exclude blocks
                                            # occur within larger legimtimate blocks. If we don't make sure that the
                                            # overlap is indicative a good match than we will lose a lot of good data.

    new_bed = ["\t".join([l[0],l[1],l[2]]) for l in bed_obj]
    return new_bed

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










def parse_regions(in_file,out_stem,keep_bed,exclude_bed):
    f = open(in_file)
    lines = f.read()
    lines = lines.rstrip("\n").split("\n")
    raw_beds = sep_valid_beds(lines)
    raw_beds['A'] = filter_bed(raw_beds['A'],keep_bed,exclude_bed)
    raw_beds['B'] = filter_bed(raw_beds['B'], keep_bed, exclude_bed)
    beds = create_merged_bedobject(raw_beds,600000)
    beds, distance = check_for_overlap(raw_beds, beds, 600000)
    sort_beds(beds,out_stem)





if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="""
    This script takes a series of region calls and filters to include targeted bed regions or exclude targeted bed 
    regions. 
    """)
    parser.add_argument('-i','--infile',
                        help='region_calls file produced by ab_call.py',
                        required=True,
                        type=str,
                        dest='in_file')
    parser.add_argument('-o','--outfile',
                        required=True,
                        type=str,
                        dest="out")
    parser.add_argument('-k','--regions_to_keep',
                        required=False,
                        default=None,
                        type=str,
                        dest='keep')
    parser.add_argument('-e', '--regions_to_exclude',
                        required=False,
                        default=None,
                        type=str,
                        dest='exclude')
    argv = parser.parse_args()



    parse_regions(argv.in_file,
                  argv.out,
                  argv.keep,
                  argv.exclude)
