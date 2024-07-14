import re
from uuid import uuid4
from collections import defaultdict
from typing import List
from .models import GenomeRange, GenomePosition

class Gene:
    def __init__(self,name,gene_id,strand,chrom,start,end,length):
        self.name = name
        self.gene_id = gene_id
        self.strand = strand
        self.chrom = chrom
        self.feature_start = start
        self.feature_end = end
        self.start = self.feature_start if strand=="+" else self.feature_end
        self.end = self.feature_end if strand=="+" else self.feature_start
        self.length = length
        self.transcripts = []
    def __repr__(self) -> str:
        return f"Gene: {vars(self)}"
    def __contains__(self,position: tuple[str,int]) -> bool:
        chrom,pos = position
        pos = GenomePosition(chrom=chrom,pos=pos)
        range = GenomeRange(chrom=self.chrom,start=self.feature_start,end=self.feature_end)
        return pos in range

class Exon:
    def __init__(self, chrom, start, end, strand, phase):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.phase = phase
    def __repr__(self):
        return "Exon: %s-%s (%s)" % (self.start, self.end, self.strand)

class Transcript:
    def __init__(self,name):
        self.name = name
        self.exons = []

def load_gff(gff) -> List[Gene]:
    GFF = open(gff)
    genes = {}
    relationships = {}
    items = {}
    while True:
        l = GFF.readline().strip()
        if not l: break
        if l[0]=="#": continue
        if l.strip()=='': continue
        fields = l.rstrip().split("\t")
        strand = fields[6]
        chrom = fields[0]
        p1 = int(fields[3])
        p2 = int(fields[4])
        feature_id = re.search("ID=([^;]*)",l)
        feature_id = feature_id.group(1) if feature_id else None
        parent_id = re.search("Parent=([^;]*)",l)
        parent_id = parent_id.group(1) if parent_id else None
        
        if fields[2] in ["gene","pseudogene","rRNA_gene","ncRNA_gene","protein_coding_gene","tRNA_gene"]:
            gene_length = p2-p1+1
            
            gene_id = None
            search_strings = [
                "ID=gene:([a-zA-Z0-9\.\-\_]+)",
                "gene_id=([a-zA-Z0-9\.\-\_]+)",
                "ID=([a-zA-Z0-9\.\-\_]+)",
                "locus_tag=([a-zA-Z0-9\.\-\_]+)",
            ]
            for s in search_strings:
                re_obj = re.search(s,l)
                if re_obj:
                    gene_id = re_obj.group(1)
                    break
            if not gene_id:
                continue
            re_obj = re.search("Name=([a-zA-Z0-9\.\-\_\(\)]+)",l)
            gene_name = re_obj.group(1) if re_obj else gene_id
            start = p1
            end =  p2
            
            items[feature_id] = Gene(gene_name,gene_id,strand,chrom,start,end,gene_length)
            relationships[feature_id] = None

        if fields[2] in ["mRNA","transcript"]:
            items[feature_id] = Transcript(feature_id)
            relationships[feature_id] = parent_id

        if fields[2] in ["CDS"]:
            if fields[7]=="":
                continue
            phase = int(fields[7])
            _id = str(uuid4())
            items[_id] = Exon(chrom,p1,p2,strand,phase)
            relationships[_id] = parent_id

    transcript_exons = defaultdict(list)
    for item in items:
        if isinstance(items[item],Exon):
            transcript_exons[relationships[item]].append(items[item])

    for item in items:
        if isinstance(items[item],Transcript):
            exons = transcript_exons[item]
            items[item].exons = sorted(exons,key=lambda x: x.start)
            gene = items[relationships[item]]
            gene.transcripts.append(items[item])

    for item in items:
        if isinstance(items[item],Gene):
            genes[items[item].gene_id] = items[item]


    return list(genes.values())