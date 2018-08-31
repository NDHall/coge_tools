"""
Here I am putting together and testing
pybedtools.
"""
import pybedtools
a = pybedtools.example_bedtool('a.bed')
b = pybedtools.example_bedtool('b.bed')

for interval in a :
    print("iteraring:",a[0])
a_and_b = a.intersect(b)

c = a.cat(b,
          postmerge=False)
#print(a)
#print(b)
#print(a_and_b)
#print(b)
#print(a)
#print(c.fn)
#print(c)
# this can be used update list of regions.
# Will probably need cleaned up afterward.
bed4 = pybedtools.BedTool('chrX 1   100 someTag',from_string=True)
bed4b = pybedtools.BedTool('chrX  1   1000    someTag', from_string=True)
bed4_total = bed4.cat(bed4b,postmerge=False)

"""
wo=True preserves both lines from previous bed and allows us 
to see matches.
"""

bed_inter = bed4.intersect(bed4b, wo=True )

"""
This can be used to determine if coordinates overlap with set.
"""
print(str(bed_inter[0]).split("\t")[4:])
print(len(bed4_total))
"""
Below statment can be used to create unique IDs
"""
print("_".join(str(bed4[0]).split("\t")))
#print(bed4_total)
