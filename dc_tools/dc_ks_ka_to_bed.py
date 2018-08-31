#ToDo: create new ks/ka syn class and region class
import dc_parse
class Ks_syn(dc_parse.Syn):

    def __init__(self,block_num, score, org_a, org_b, orient, num_gene_pairs,ks,ka ):
        super().__init__(block_num, score, org_a, org_b, orient, num_gene_pairs)
        self.ks = float(ks)
        self.ka = float(ka)

#ToDo: create ks/ka parsing function.
#ToDo: make ks/ka bed producer.
#ToDo: ks/ka filter
