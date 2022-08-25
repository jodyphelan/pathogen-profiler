from tqdm import tqdm 
import statistics as stats
from .utils import debug, infolog,revcom
from itertools import combinations, product
import os
def get_canonical_kmer(kmer):
    t = {
        "A":1,
        "C":2,
        "G":3,
        "T":4
    }
    rkmer = revcom(kmer)
    nkmer = int("".join([str(t[x]) for x in list(kmer)]))
    nrkmer = int("".join([str(t[x]) for x in list(rkmer)]))
    return kmer if nkmer<nrkmer else rkmer

def mutate_kmer(kmer,d=1):

    def generate(s, d=1):
        N = len(s)
        letters = 'ACGT'
        pool = list(s)

        for indices in combinations(range(N), d):
            for replacements in product(letters, repeat=d):
                skip = False
                for i, a in zip(indices, replacements):
                    if pool[i] == a: skip = True
                if skip: continue

                keys = dict(zip(indices, replacements))
                yield ''.join([pool[i] if i not in indices else keys[i] 
                            for i in range(N)])

    kmers = set([get_canonical_kmer(k) for k in [get_canonical_kmer(kmer)] + list(generate(kmer,d=d))])
    return kmers

class kmer_dump:
    def __init__(self,kmer_file):
        self.kmer_file = kmer_file


    def load_kmer_counts(self,kmer_db_file,remove_after_processing=True):
        self.kmer_counts = []
        kmers = {}
        for l in open(kmer_db_file):
            row = l.strip().split("\t")
            for k in mutate_kmer(row[0],d=1):
                kmers[k] = row[1]
        tmp_counts = {}
        infolog(f"Looking for {len(kmers)} kmers")
        for l in tqdm(open(self.kmer_file)):
            row = l.strip().split()
            if row[0] not in kmers: continue
            tmp_counts[row[0]] = int(row[1])
        
        for l in open(kmer_db_file):
            row = l.strip().split("\t")
            count = sum([tmp_counts.get(k,0) for k in mutate_kmer(row[0],d=1)])
            self.kmer_counts.append({"name":kmers[get_canonical_kmer(row[0])],"seq":row[0],"count":count})
        if remove_after_processing:
            os.remove(self.kmer_file)
        return self.kmer_counts

    def get_taxonomic_support(self,kmer_db_file,output_kmer_counts=None):
        if not hasattr(self, 'kmer_counts'):
            self.load_kmer_counts(kmer_db_file)
        
        if output_kmer_counts:
            with open(output_kmer_counts,"w") as O:
                for k in self.kmer_counts:
                    O.write("%(name)s\t%(seq)s\t%(count)s\n" % k)
        
        
        taxon_set = set(l.strip().split("\t")[1] for l in open(kmer_db_file))
        
        taxon_support = []
        for s in taxon_set:
            support = [x["count"] for x in self.kmer_counts if x["name"]==s]
            if len([x for x in support if x!=0])<len(support)/2:
                continue
            mean = stats.mean(support)
            std = stats.stdev(support)
            taxon_support.append({"species":s,"mean":mean,"std":std})

        return taxon_support
