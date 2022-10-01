import csv
from glob import glob
import json
import re
from collections import defaultdict
import sys
from datetime import datetime
from .utils import load_gff, run_cmd, cmd_out, errlog, unlist, debug, infolog, successlog
from .fasta import fasta
import os
import shutil
from uuid import uuid4
import pathogenprofiler as pp


supported_so_terms = [
    'coding_sequence_variant', 'chromosome', 'duplication', 'inversion', 'coding_sequence_variant', 
    'inframe_insertion', 'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion', 
    'downstream_gene_variant', 'exon_variant', 'exon_loss_variant', 'exon_loss_variant', 'duplication', 
    'duplication', 'inversion', 'inversion', 'frameshift_variant', 'gene_variant', 'feature_ablation', 
    'duplication', 'gene_fusion', 'gene_fusion', 'bidirectional_gene_fusion', 'rearranged_at_DNA_level', 
    'intergenic_region', 'conserved_intergenic_variant', 'intragenic_variant', 'intron_variant', 
    'conserved_intron_variant', 'miRNA', 'missense_variant', 'initiator_codon_variant', 'stop_retained_variant', 
    'protein_protein_contact', 'structural_interaction_variant', 'rare_amino_acid_variant', 
    'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant', 'splice_region_variant', 
    'splice_region_variant', 'stop_lost', '5_prime_UTR_premature_', 'start_codon_gain_variant', 
    'start_lost', 'stop_gained', 'synonymous_variant', 'start_retained', 'stop_retained_variant', 
    'transcript_variant', 'transcript_ablation', 'regulatory_region_variant', 'upstream_gene_variant', 
    '3_prime_UTR_variant', '3_prime_UTR_truncation + exon_loss', '5_prime_UTR_variant', 
    '5_prime_UTR_truncation + exon_loss_variant', 'sequence_feature + exon_loss_variant'
]

def generate_kmer_database(kmer_file,outfile):
    from itertools import combinations, product

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

    with open(outfile,"w") as O:
        for l in open(kmer_file):
            row = l.strip().split()
            kmers = [row[0]] + list(generate(row[0]))
            O.write("%s\t%s\t%s\n" % (row[1],len(kmers),"\t".join(kmers)))



def fa2dict(filename):
    fa_dict = {}
    seq_name = ""
    for l in open(filename):
        line = l.rstrip()
        if line=="":continue
        if line[0] == ">":
            seq_name = line[1:].split()[0]
            fa_dict[seq_name] = []
        else:
            fa_dict[seq_name].append(line)
    result = {}
    for seq in fa_dict:
        result[seq] = "".join(fa_dict[seq])
    return result

def revcom(s):
    """Return reverse complement of a sequence"""
    def complement(s):
                    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                    letters = list(s)
                    letters = [basecomplement[base] for base in letters]
                    return ''.join(letters)
    return complement(s[::-1])


def extract_genome_positions(db,gene):
    pos = []
    for mut in db[gene]:
        if any([a["type"]=="drug" for a in db[gene][mut]["annotations"]]):
            if mut[:1] not in ["p.","c.","n."]: continue
            pos.extend(db[gene][mut]["genome_positions"])
    return list(set(pos))

def write_bed(db,gene_dict,gene_info,ref_fasta,outfile,padding=200):

    lines = []
    for gene in gene_dict:
        if gene not in gene_info:
            errlog("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
            quit()
        if gene_info[gene].locus_tag in db:
            genome_positions = extract_genome_positions(db,gene_info[gene].locus_tag)
            if len(genome_positions)>0 and (gene_info[gene].feature_start > min(genome_positions)):
                genome_start = min(genome_positions) - padding
            else:
                genome_start = gene_info[gene].feature_start - padding
            
            if len(genome_positions)>0 and (gene_info[gene].feature_end < max(genome_positions)):
                genome_end = max(genome_positions) + padding
            else:
                genome_end = gene_info[gene].feature_end + padding
        else:
            genome_start = gene_info[gene].feature_start - padding
            genome_end = gene_info[gene].feature_end + padding

        if genome_start<1:
            genome_start = 1
        if genome_end>len(ref_fasta[gene_info[gene].chrom]):
            genome_end = len(ref_fasta[gene_info[gene].chrom])

        drugs = [d for d in gene_dict[gene] if d!=""]
        if len(drugs)==0:
            drugs = "None"
        else:
            drugs = ",".join(drugs)
        lines.append([
            gene_info[gene].chrom,
            str(genome_start),
            str(genome_end),
            gene_info[gene].locus_tag,
            gene_info[gene].name,
            drugs
        ])
    with open(outfile,"w") as O:
        for line in sorted(lines,key=lambda x: (x[0],int(x[1]))):
            O.write("%s\n" %"\t".join(line))

# def load_gene_info(filename):
#     gene_info = {}
#     for l in open(filename):
#         row = l.rstrip().split()
#         strand = "-" if row[0][-1]=="c" else "+"
#         gene_info[row[0]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
#         gene_info[row[1]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
#     return gene_info

def get_ann(variants,snpEffDB):
    uuid = str(uuid4()) #"463545ef-71fc-449b-8f4e-9c907ee6fbf5"
    with open(uuid,"w") as O:
        O.write('##fileformat=VCFv4.2\n')
        O.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        O.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest\n')
        for var in variants.values():
            O.write("%(chrom)s\t%(pos)s\t.\t%(ref)s\t%(alt)s\t255\t.\t.\tGT\t1\n" % var) 
    results = {}
    keys = list(variants.keys())
    vals = list(variants.values())
    i = 0
    for l in cmd_out(f"snpEff ann {snpEffDB} {uuid}"):
        if l[0]=="#": continue
        row = l.strip().split()
        for ann in row[7].split(","):
            a = ann.split("|")
            if len(a)!=16:continue

            if vals[i]["gene"] in [a[3],a[4]]:
                results[keys[i]] = a[9] if vals[i]["type"]=="nucleotide" else a[10]
        i+=1
    os.remove(uuid)
    return results

def assign_gene_to_amplicon(genes,chrom,start,end):
    l = []
    for g in genes.values():
        if g.chrom!=chrom: continue
        overlap = set(range(g.feature_start,g.feature_end)).intersection(set(range(int(start),int(end))))
        if overlap:
            l.append((g.locus_tag,g.name,len(overlap)))
    return tuple(sorted(l,key=lambda x:x[2],reverse=True)[0][:2])



def assign_amplicon_drugs(db,chrom,start,end):
    d = set()
    for gene in db:
        for change in db[gene]:
            if db[gene][change]['chromosome']!=chrom: continue
            if db[gene][change]['genome_positions']==None: continue
            if set(db[gene][change]['genome_positions']).intersection(set(range(int(start),int(end)))):
                for ann in db[gene][change]['annotations']:
                    if "drug" in ann:
                        d.add(ann['drug'])
    return d

def write_amplicon_bed(ref_seq,genes,db,primer_file,outfile):
    ref = fasta(ref_seq)
    with open(outfile,"w") as O:
        for chrom,start,end,amplicon_name in ref.get_amplicons(primer_file):
            locus_tag,gene_name = assign_gene_to_amplicon(genes,chrom,start,end)
            drugs = ",".join(assign_amplicon_drugs(db,chrom,start,end))
            if drugs=="":
                drugs = "None"
            O.write(f"{chrom}\t{start}\t{end}\t{locus_tag}\t{gene_name}\t{drugs}\t{amplicon_name}\n")

def get_snpeff_formated_mutation_list(csv_file,ref,gff,snpEffDB):
    genes = load_gff(gff,aslist=True)
    refseq = fa2dict(ref)
    mutations  =  {}
    converted_mutations = {}
    for row in csv.DictReader(open(csv_file)):
        gene = [g for g in genes if g.name==row["Gene"] or g.locus_tag==row["Gene"]][0]
        r = re.search("n.([0-9]+)([ACGT]+)>([ACGT]+)",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
            
        r = re.search("n.([0-9]+)([ACGT]+)>([ACGT]+)",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = f"n.{r.group(1)}{r.group(2).upper()}>{r.group(3).upper()}"
        r = re.search("p\..+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        r = re.search("c.-[0-9]+[ACGT]>[ACGT]",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]

        r = re.search("c.[0-9]+dup[ACGT]+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        r = re.search("c.[0-9]+_[0-9]+dup[ACGT]+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        



        r = re.search("c.([0-9]+)del",row["Mutation"])
        if r:
            # "ethA" "c.341del"
            del_start = int(r.group(1))
            del_end = int(r.group(1))
            if gene.strand == "+":
                # rpoB "c.1282_1290del"
                genome_start = gene.start + del_start - 2
                genome_end = gene.start + del_end 
            else:
                # "ethA" "c.1057_1059del"
                genome_start = gene.start - del_end
                genome_end = gene.start - del_start + 2
            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        r = re.search("c.([0-9]+)_([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
                # rpoB "c.1282_1290del"
                genome_start = gene.start + del_start - 2
                genome_end = gene.start + del_end 
            else:
                # "ethA" "c.1057_1059del"
                genome_start = gene.start - del_end
                genome_end = gene.start - del_start + 2
            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        r = re.search("c.-([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(1))
            if gene.strand == "+":
               # "embA" "c.-29_-28del"
                genome_start = gene.start - del_start - 1
                genome_end = gene.start - del_end + 1
            else:
                # "alr" "c.-283_-280delCAAT"
                genome_start = gene.start + del_end - 1
                genome_end = gene.start + del_start + 1

            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        
        r = re.search("c.(-[0-9]+)_(-[0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
               # "embA" "c.-29_-28del"
                genome_start = gene.start + del_start - 1
                genome_end = gene.start + del_end + 1
            else:
                # "alr" "c.-283_-280delCAAT"
                genome_start = gene.start - del_end - 1
                genome_end = gene.start - del_start + 1
            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        
        r = re.search("c.(-[0-9]+)_([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
               # "ethA" "c.-1058_968del"
                genome_start = gene.start + del_start -1
                genome_end = gene.start + del_end 
                # quit("Need to define!")

            else:
               # "ethA" "c.-1058_968del"
                genome_start = gene.start - del_end 
                genome_end = gene.start - del_start + 1

            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}


        r = re.search("c.([0-9]+)_([0-9]+)ins([ACGT]+)", row["Mutation"])
        if r:
            ins_start = int(r.group(1))
            ins_end = int(r.group(2))
            ins_seq = r.group(3)
            if gene.strand == "+":
                # "rpoB" "c.1296_1297insTTC"
                genome_start = gene.start + ins_start - 1 
                genome_end = gene.start + ins_end - 1
            else:
                # "pncA" "c.521_522insT"
                ins_seq = pp.revcom(ins_seq)
                genome_start = gene.start - ins_start 
                genome_end = gene.start - ins_end + 2

            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref + ins_seq
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}
        

        r = re.search("c.(-[0-9]+)_(-[0-9]+)ins([ACGT]+)",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            ins_seq = r.group(3)
            if gene.strand == "+":
               # "rrs" "c.-29_-28insATAC"
                genome_start = gene.start + del_start 
                genome_end = gene.start + del_end 
            else:
                # "alr" "c.-283_-280delCAAT"
                ins_seq = pp.revcom(ins_seq)
                genome_start = gene.start - del_end 
                genome_end = gene.start - del_start 
            ref = refseq[gene.chrom][genome_start-1:genome_end-1]
            alt = ref + ins_seq
            mutations[(row["Gene"],row["Mutation"])] = {"chrom":gene.chrom,"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}


        if (row["Gene"],row["Mutation"]) not in converted_mutations and (row["Gene"],row["Mutation"]) not in mutations:
                if row['Mutation'] in supported_so_terms:
                    converted_mutations[(row["Gene"],row['Mutation'])] = row['Mutation']
        if (row["Gene"],row["Mutation"]) not in converted_mutations and (row["Gene"],row["Mutation"]) not in mutations:
            raise Exception(f"Don't know how to handle this mutation: {row['Gene']} {row['Mutation']}")            

    infolog("Converting %s mutations" % len(mutations))
    if len(mutations)>0:
        mutation_conversion = get_ann(mutations,snpEffDB)
        for key in mutation_conversion:
            converted_mutations[key] = mutation_conversion[key]
    return converted_mutations
    



def get_genome_position(gene_object,change):
    for term in supported_so_terms:
        if term in change:
            return None
    if "any_missense_codon" in change:
        codon = int(change.replace("any_missense_codon_",""))
        change = f"p.Xyz{codon}Xyz"


    g = gene_object
    r = re.search("p.[A-Za-z]+([0-9]+)",change)
    if r:
        codon = int(r.group(1))
        if g.strand=="+":
            p = g.start + codon*3-3
            return [p,p+1,p+2]
        else:
            p = g.start - codon*3 + 1
            return [p,p+1,p+2]
    r = re.search("c.(-[0-9]+)[ACGT]+>[ACGT]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos
            return [p]
        else:
            p = g.start - pos
            return [p]
    r = re.search("n.([0-9]+)[ACGT]+>[ACGT]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos -1
            return [p]
    r = re.search("[nc].([\-\*0-9]+)_([\-\*0-9]+)ins[A-Z]+",change)
    if r:
        if "*" in r.group(1):
            if g.strand=="+":
                pos = g.feature_end - g.feature_start + int(r.group(1).replace("*","")) + 1
            else:
                pos = g.feature_end - g.feature_start - int(r.group(1).replace("*","")) - 1
        else:
            pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos -1
            return [p, p+1]
        else:
            p = g.start - pos 
            return [p, p+1]

    r = re.search("[nc].([\-\*0-9]+)_([\-\*0-9]+)del[A-Z]*",change)
    if r:
        if "*" in r.group(1):
            if g.strand=="+":
                pos1 = g.feature_end - g.feature_start + int(r.group(1).replace("*",""))
            else:
                pos1 = g.feature_end - g.feature_start - int(r.group(1).replace("*",""))
        else:
            pos1 = int(r.group(1))
        if "*" in r.group(2):
            if g.strand=="+":
                pos2 = g.feature_end - g.feature_start + int(r.group(2).replace("*",""))
            else:
                pos2 = g.feature_end - g.feature_start - int(r.group(2).replace("*",""))
        else:
            pos2 = int(r.group(2))
        if g.strand=="+":
            p1 = g.start + pos1 -1
            p2 = g.start + pos2 -1
            if pos1<0:
                p1+=1
                p2+=1
            return list(range(p1,p2+1))
        else:
            p1 = g.start - pos1 + 1
            p2 = g.start - pos2 + 1
            if pos1<0:
                p1-=1
            if pos2<0:
                p2-=1
            return list(range(p2,p1+1))

    r = re.search("[nc].([\-\*0-9]+)del[A-Z]+",change)
    if r:
        if "*" in r.group(1):
            if g.strand=="+":
                pos = g.feature_end - g.feature_start + int(r.group(1).replace("*","")) + 1
            else:
                pos = g.feature_end - g.feature_start - int(r.group(1).replace("*","")) - 1
        else:
            pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]

    r = re.search("[nc].([\-0-9]+)dup[A-Z]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]
    
    r = re.search("[nc].([\-0-9]+)_([\-0-9]+)dup[A-Z]+",change)
    if r:
        pos1 = int(r.group(1))
        pos2 = int(r.group(2))
        if g.strand=="+":
            p1 = g.start + pos1 - 1
            p2 = g.start + pos2 - 1
            return list(range(p1,p2+1))
        else:
            p1 = g.start - pos1 + 1
            p2 = g.start - pos2 + 1
            return list(range(p2,p1+1))

    

    # r = re.search("n.([0-9]+)([0-9]+)dup[A-Z]+",change)
    # if r:
    #     pos1 = int(r.group(1))
    #     pos2 = int(r.group(2))
    #     if g.strand=="+":
    #         p1 = g.start + pos1 - 1
    #         p2 = g.start + pos2 - 1
    #         return list(range(p1,p2+1))
    #     else:
    #         p = g.start - pos + 1
    #         quit(f"Don't know how to handle {change}")
    #         return [p]
    quit(f"Don't know how to handle {str(vars(g))} {change}")

def match_ref_chrom_names(source,target):
    source_fa = fa2dict(source)
    source_fa_size = {s:len(source_fa[s]) for s in source_fa}
    target_fa = fa2dict(target)
    target_fa_size = {s:len(target_fa[s]) for s in target_fa}
    conversion = {}
    for s in target_fa:
        tlen = target_fa_size[s]
        tmp = [x[0] for x in source_fa_size.items() if x[1]==tlen]
        if len(tmp)==1:
            conversion[s] = tmp[0]
    return conversion


def create_db(args,extra_files = None):
    variables = json.load(open("variables.json"))    
    genome_file = "%s.fasta" % args.prefix
    gff_file = "%s.gff" % args.prefix
    bed_file = "%s.bed" % args.prefix
    json_file = "%s.dr.json" % args.prefix
    version_file = "%s.version.json" % args.prefix

    if not extra_files:
        extra_files = {}

    if args.match_ref:
        chrom_conversion = match_ref_chrom_names(args.match_ref,"genome.fasta")
        shutil.copyfile(args.match_ref,genome_file)
    else:
        chrom_conversion = match_ref_chrom_names("genome.fasta","genome.fasta")
        shutil.copyfile("genome.fasta",genome_file)    
    
    with open(gff_file,"w") as O:
        for l in open("genome.gff"):
            if l.strip()=="": continue
            if l[0]=="#":
                O.write(l)
            else:
                row = l.strip().split()
                if row[0] in chrom_conversion:
                    row[0] = chrom_conversion[row[0]]
                    O.write("\t".join(row)+"\n")        

    genes = load_gff(gff_file)
    gene_name2gene_id = {g.name:g.locus_tag for g in genes.values()}
    gene_name2gene_id.update({g.locus_tag:g.locus_tag for g in genes.values()})
    db = {}
    locus_tag_to_drug_dict = defaultdict(set)
    with open(args.prefix+".conversion.log","w") as L:
        if args.csv:
            mutation_lookup = get_snpeff_formated_mutation_list(args.csv,"genome.fasta","genome.gff",json.load(open("variables.json"))["snpEff_db"])
            for row in csv.DictReader(open(args.csv)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                drug = row["Drug"].lower()
                mut = mutation_lookup[(row["Gene"],row["Mutation"])]
                if args.include_original_mutation:
                    row["original_mutation"] = row["Mutation"]
                if mut!=row["Mutation"]:
                    L.write(f"Converted {row['Gene']} {row['Mutation']} to {mut}\n")
                locus_tag_to_drug_dict[locus_tag].add(drug)
                if locus_tag not in db:
                    db[locus_tag] = {}
                if mut not in db[locus_tag]:
                    db[locus_tag][mut] = {"annotations":[]}

                tmp_annotation = {"type":"drug","drug":row["Drug"]}
                annotation_columns = set(row.keys()) - set(["Gene","Mutation","Drug"])
                for col in annotation_columns:
                    if row[col]=="":continue
                    tmp_annotation[col.lower()] = row[col]
                db[locus_tag][mut]["annotations"].append(tmp_annotation)
                db[locus_tag][mut]["genome_positions"] = get_genome_position(genes[locus_tag],mut) if mut not in supported_so_terms else None
                db[locus_tag][mut]["chromosome"] = genes[locus_tag].chrom
        if args.other_annotations:
            mutation_lookup = get_snpeff_formated_mutation_list(args.other_annotations,"genome.fasta","genome.gff",json.load(open("variables.json"))["snpEff_db"])
            for row in csv.DictReader(open(args.other_annotations)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                mut = mutation_lookup[(row["Gene"],row["Mutation"])]
                if mut!=row["Mutation"]:
                    L.write(f"Converted {row['Gene']} {row['Mutation']} to {mut}\n")
                if locus_tag not in db:
                    db[locus_tag] = {}
                if mut not in db[locus_tag]:
                    db[locus_tag][mut] = {"annotations":[]}
                tmp_annotation = {"type":row["Type"]}
                if args.include_original_mutation:
                    tmp_annotation["original_mutation"] = row["Mutation"]


                for x in row["Info"].split(";"):
                    key,val = x.split("=")
                    tmp_annotation[key.lower()] = val
                    if key=="drug":
                        locus_tag_to_drug_dict[locus_tag].add(val)
                db[locus_tag][mut]["annotations"].append(tmp_annotation)
                db[locus_tag][mut]["genome_positions"] = get_genome_position(genes[locus_tag],mut)
                db[locus_tag][mut]["chromosome"] = genes[locus_tag].chrom

        if args.watchlist:
            for row in csv.DictReader(open(args.watchlist)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                for d in row["Drug"].split(","):
                    drug = d.lower()
                    locus_tag_to_drug_dict[locus_tag].add(drug)


        version = {"name":args.prefix}
        if not args.custom:
            for l in cmd_out("git log | head -4"):
                row = l.strip().split()
                if row == []: continue
                version[row[0].replace(":","")] = " ".join(row[1:])
            version["commit"] = version["commit"][:7]
        else:
            version["Date"] = str(datetime.now()) if not args.db_date else args.db_date
            version["name"] = args.db_name if args.db_name else "NA"
            version["commit"] = args.db_commit if args.db_name else "NA"
            version["Author"] = args.db_author if args.db_author else "NA"

        json.dump(version,open(version_file,"w"))
        json.dump(db,open(json_file,"w"))
        
        
        for file in extra_files.values():
            target = f"{args.prefix}.{file}"
            shutil.copyfile(file,target)

        
        if "barcode" in extra_files:
            barcode_file = f"{args.prefix}.{extra_files['barcode']}"

            with open(barcode_file,"w") as O:
                for l in open("barcode.bed"):
                    if l[0]=="#": continue
                    row = l.strip().split("\t")
                    row[0] = chrom_conversion[row[0]]
                    O.write("\t".join(row)+"\n")
        
        if "amplicon_primers" in vars(args) and args.amplicon_primers:
            write_amplicon_bed(genome_file,genes,db,args.amplicon_primers,bed_file)
            variables['amplicon'] = True
        else:
            ref_fasta_dict = fa2dict(genome_file)
            write_bed(db,locus_tag_to_drug_dict,genes,ref_fasta_dict,bed_file)
            variables['amplicon'] = False
        
                
        if list(chrom_conversion.keys())!=list(chrom_conversion.values()):
            variables["chromosome_conversion"] = {"target":list(chrom_conversion.keys()),"source":list(chrom_conversion.values())}
        variables_file = args.prefix+".variables.json"
        variables["files"] = {
            "ref": genome_file,
            "gff": gff_file,
            "bed": bed_file,
            "version": version_file,
            "json_db": json_file,
            "variables": variables_file
        }
        if extra_files:
            for key,val in extra_files.items():
                variables["files"][key] = f"{args.prefix}.{val}"
        json.dump(variables,open(variables_file,"w"))
        
        if os.path.isfile("snpEffectPredictor.bin"):
            snpeff_db_name = json.load(open("variables.json"))["snpEff_db"]
            load_snpEff_db("snpEffectPredictor.bin",snpeff_db_name)
        
        if args.load:
            load_db(variables_file,args.software_name)

def load_db(variables_file,software_name,source_dir="."):
    variables = json.load(open(variables_file))
    load_dir = f"{sys.base_prefix}/share/{software_name}"
    if not os.path.isdir(load_dir):   
        os.mkdir(load_dir)

    for key,val in variables['files'].items():
        source = f"{source_dir}/{val}"
        target = f"{load_dir}/{val}"
        infolog(f"Copying file: {source} ---> {target}")
        shutil.copyfile(source,target)
        if key=="ref":
            pp.run_cmd(f"bwa index {target}")
            pp.run_cmd(f"samtools faidx {target}")
            tmp = target.replace(".fasta","")
            pp.run_cmd(f"samtools dict {target} -o {tmp}.dict")
    
    successlog("Sucessfully imported library")

def get_variable_file_name(software_name,library_name):
    library_prefix = f"{sys.base_prefix}/share/{software_name}/{library_name}"
    return f"{library_prefix}.variables.json"

def get_db(software_name,db_name):
    if "/" in db_name:
        share_path = "/".join(db_name.split("/")[:-1])
        db_name = db_name.split("/")[-1]
        variable_file_name = f"{share_path}/{db_name}.variables.json"
    else:
        share_path = f"{sys.base_prefix}/share/{software_name}/"
        variable_file_name = get_variable_file_name(software_name,db_name)
    
    if not os.path.isfile(variable_file_name):
        return None
    variables = json.load(open(variable_file_name))
    for key,val in variables['files'].items():
        infolog(f"Using {key} file: {share_path}/{val}")
        if ".json" in val:
            variables[key] = json.load(open(f"{share_path}/{val}"))
        else:
            variables[key] = f"{share_path}/{val}"
    
    return variables    

def list_db(software_name):
    share_path = f"{sys.base_prefix}/share/{software_name}"
    if not os.path.isdir(share_path):
        return []
    return [json.load(open(f"{share_path}/{f}")) for f in os.listdir(share_path) if f.endswith(".version.json")]



def create_species_db(args,extra_files = None):
    if not extra_files:
        extra_files = {}
    version = {"name":args.prefix}
    if not args.db_name:
        for l in pp.cmd_out("git log | head -4"):
            row = l.strip().split()
            if row == []: continue
            version[row[0].replace(":","")] = " ".join(row[1:])
        version["commit"] = version["commit"][:7]
    else:
        version["Date"] = str(datetime.now()) if not args.db_date else args.db_date
        version["name"] = args.db_name if args.db_name else "NA"
        version["commit"] = args.db_commit if args.db_name else "NA"
        version["Author"] = args.db_author if args.db_author else "NA"

    kmer_file = args.prefix+".kmers.txt"
    version_file = args.prefix+".version.json"
    shutil.copyfile(args.kmers,kmer_file)
    json.dump(version,open(version_file,"w"))
    for file in extra_files.values():
            target = f"{args.prefix}.{file}"
            shutil.copyfile(file,target)
    variables_file = args.prefix+".variables.json"
    variables = {}
    variables["files"] = {
        "kmers": kmer_file,
        "version": version_file,
        "variables": variables_file
    }

    if extra_files:
        for key,val in extra_files.items():
            variables["files"][key] = f"{args.prefix}.{val}"
    json.dump(variables,open(variables_file,"w"))

    if args.load:
        load_dir = f"{sys.base_prefix}/share/{args.software_name}"
        if not os.path.isdir(load_dir):   
            os.mkdir(load_dir)

        for key,val in variables['files'].items():
            target = f"{load_dir}/{val}"
            infolog(f"Copying file: {val} ---> {target}")
            shutil.copyfile(val,target)

def get_snpeff_dir():
    tmp = glob(f"{sys.base_prefix}/share/*snpeff*")
    if len(tmp)>0:
        return tmp[0]
    else: 
        return None

def load_snpEff_db(bin_file,genome_name):
    snpeff_dir = get_snpeff_dir()
    snpeff_config = f"{snpeff_dir}/snpEff.config"
    with open(snpeff_config,"a") as F:
        F.write(f"\n{genome_name}.genome : {genome_name}\n")
    
    data_dir = f"{snpeff_dir}/data"
    genome_dir = f"{data_dir}/{genome_name}"
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    if not os.path.isdir(genome_dir):
        os.mkdir(genome_dir)
    shutil.copyfile(bin_file,f"{genome_dir}/{bin_file}")
