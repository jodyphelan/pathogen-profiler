from tqdm import tqdm 
import statistics as stats
from .utils import debug, infolog,revcom
from itertools import combinations, product
import os




class kmer_dump:
    def __init__(self,kmer_file,counter):
        self.kmer_file = kmer_file
        if counter=="kmc":
            nuc_order = "ACGT"
        elif counter=="dsk":
            nuc_order = "ACTG"
        self.nuc_order = {n:i for i,n in enumerate(nuc_order)}

    def load_kmer_counts(self,kmer_db_file,remove_after_processing=True,max_mismatch=1):
        self.kmer_counts = []
        kmers = {}
        for l in open(kmer_db_file):
            row = l.strip().split("\t")
            for k in self.mutate_kmer(row[0],d=max_mismatch):
                kmers[k] = row[1]
        tmp_counts = {}
        infolog(f"Looking for {len(kmers)} kmers")
        for l in tqdm(open(self.kmer_file)):
            row = l.strip().split()
            if row[0] not in kmers: continue
            tmp_counts[row[0]] = int(row[1])
        
        for l in open(kmer_db_file):
            row = l.strip().split("\t")
            count = sum([tmp_counts.get(k,0) for k in self.mutate_kmer(row[0],d=max_mismatch)])
            self.kmer_counts.append({"name":kmers[self.get_canonical_kmer(row[0])],"seq":row[0],"count":count})
        if remove_after_processing:
            os.remove(self.kmer_file)
        return self.kmer_counts

    def get_taxonomic_support(self,kmer_db_file,output_kmer_counts=None):
        if not hasattr(self, 'kmer_counts'):
            self.load_kmer_counts(kmer_db_file,max_mismatch=0)
        
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

    def get_canonical_kmer(self,kmer):
        t = self.nuc_order
        rkmer = revcom(kmer)
        nkmer = int("".join([str(t[x]) for x in list(kmer)]))
        nrkmer = int("".join([str(t[x]) for x in list(rkmer)]))
        return kmer if nkmer<nrkmer else rkmer
    
    def mutate_kmer(self,kmer,d=1):

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

        kmers = set([self.get_canonical_kmer(k) for k in [self.get_canonical_kmer(kmer)] + list(generate(kmer,d=d))])
        return kmers