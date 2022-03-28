from tqdm import tqdm 
import statistics as stats
from .utils import debug,revcom

def get_canonical_kmer(kmer):
    t = {"A":"1","C":"2","T":"3","G":"4"}
    rkmer = revcom(kmer)
    nkmer = int("".join([t[x] for x in list(kmer)]))
    nrkmer = int("".join([t[x] for x in list(rkmer)]))
    return kmer if nkmer<nrkmer else rkmer


class kmer_dump:
    def __init__(self,kmer_file):
        self.kmer_file = kmer_file


    def load_kmer_counts(self,kmer_db_file):
        self.kmer_counts = []
        kmers = {}
        for l in open(kmer_db_file):
            row = l.strip().split("\t")
            kmers[get_canonical_kmer(row[0])] = row[1]
        tmp_counts = {}
        for l in tqdm(open(self.kmer_file)):
            row = l.strip().split()
            if get_canonical_kmer(row[0]) not in kmers: continue
            tmp_counts[get_canonical_kmer(row[0])] = int(row[1])
        
        for l in open(kmer_db_file):
            row = l.strip().split("\t")
            count = tmp_counts.get(get_canonical_kmer(row[0]),0)
            self.kmer_counts.append({"name":kmers[get_canonical_kmer(row[0])],"seq":row[0],"count":count})
        return self.kmer_counts

    def get_taxonomic_support(self,kmer_db_file):
        if not hasattr(self, 'kmer_counts'):
            self.load_kmer_counts(kmer_db_file)
        
        tmp_counts = {x["seq"]:x["count"] for x in self.kmer_counts}
        
        taxon_set = set(l.strip().split("\t")[1] for l in open(kmer_db_file))
        
        taxon_support = []
        for s in taxon_set:
            support = [x["count"] for x in self.kmer_counts if x["name"]==s]
            print(s,support)
            if len([x for x in support if x!=0])<len(support)/2:
                continue
            mean = stats.mean(support)
            std = stats.stdev(support)
            taxon_support.append({"taxon":s,"mean":mean,"std":std})

        return taxon_support
