import csv
from glob import glob
import json
import re
from collections import defaultdict
import sys
from .gff import load_gff, Gene
from .fasta import Fasta
import os
import shutil
import pathogenprofiler as pp
import math
import logging
from .hgvs import verify_mutation_list
from pysam import FastaFile
from typing import List
import argparse
import subprocess as sp


def get_git_repo_info():
    git_info = {
        'repo': None,
        'branch': None,
        'commit': None,
        'status': None,
        'author': None,
        'date': None
    }
    status = sp.run('git status', shell=True, capture_output=True).stderr.decode("utf-8")

    if status.startswith('fatal'):
        log = None
        logging.warning("Not a git repository")
    else:
        log = sp.run('git log -1', shell=True, capture_output=True).stdout.decode("utf-8")
        # git config --get remote.origin.url
        git_info['repo'] = sp.run('git config --get remote.origin.url', shell=True, capture_output=True).stdout.decode("utf-8").strip()
        # git rev-parse --abbrev-ref HEAD
        git_info['branch'] = sp.run('git rev-parse --abbrev-ref HEAD', shell=True, capture_output=True).stdout.decode("utf-8").strip()
        # git rev-parse --short=8 HEAD
        git_info['commit'] = sp.run('git rev-parse --short=8 HEAD', shell=True, capture_output=True).stdout.decode("utf-8").strip()
        # git status
        tmp = sp.run('git status', shell=True, capture_output=True).stdout.decode("utf-8").strip()
        if 'modified' in tmp:
            git_info['status'] = 'modified'
        else:
            git_info['status'] = 'clean'
        # git log -1 --pretty=format:"%an"
        git_info['author'] = sp.run('git log -1 --pretty=format:"%an"', shell=True, capture_output=True).stdout.decode("utf-8").strip()
        # git log -1 --pretty=format:"%ad"
        git_info['date'] = sp.run('git log -1 --pretty=format:"%ad"', shell=True, capture_output=True).stdout.decode("utf-8").strip()

    return git_info


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
    '5_prime_UTR_truncation + exon_loss_variant', 'sequence_feature + exon_loss_variant', 'functionally_normal',
    'conservative_inframe_deletion', 'conservative_inframe_insertion'
]

def generate_kmer_database(kmer_file: str,outfile: str) -> None:
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



def revcom(s: str) -> str:
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
        if mut[:2] not in ["p.","c.","n."]: continue
        pos.extend(db[gene][mut]["genome_positions"])
    return list(set(pos))

def write_bed(
        db: dict,
        gene_dict: dict,
        gene_info: List[Gene],
        ref_file: str,
        outfile:str,
        padding: int = 200,
        gene_coordinates: List[dict] = None
    ) -> None:
    if os.path.exists(ref_file+".fai"):
        os.remove(ref_file+".fai")
    ref = FastaFile(ref_file)
    lines = []
    for gene in gene_dict:
        if gene not in gene_info:
            logging.error("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
            quit()
        if gene_coordinates and gene in gene_coordinates:
            genome_start = gene_coordinates[gene]['start']
            genome_end = gene_coordinates[gene]['end']
        elif gene_info[gene].gene_id in db:
            genome_positions = extract_genome_positions(db,gene_info[gene].gene_id)
            if gene_info[gene].strand=="+":
                if len(genome_positions)>0 and (gene_info[gene].feature_start > min(genome_positions)):
                    genome_start = min(genome_positions) - padding
                else:
                    genome_start = gene_info[gene].start - padding
                
                if len(genome_positions)>0 and (gene_info[gene].feature_end < max(genome_positions)):
                    genome_end = max(genome_positions) #+ padding
                else:
                    genome_end = gene_info[gene].end #+ padding
            else:
                if len(genome_positions)>0 and (gene_info[gene].feature_start > min(genome_positions)):
                    genome_start = min(genome_positions) #- padding
                else:
                    genome_start = gene_info[gene].end #- padding
                
                if len(genome_positions)>0 and (gene_info[gene].feature_end < max(genome_positions)):
                    genome_end = max(genome_positions) + padding
                else:
                    genome_end = gene_info[gene].start + padding
        else:
            genome_start = gene_info[gene].feature_start - padding
            genome_end = gene_info[gene].feature_end + padding

        if genome_start<1:
            genome_start = 1
        
        chrom_lengths = dict(zip(ref.references, ref.lengths))
        if genome_end > chrom_lengths[gene_info[gene].chrom]:
            genome_end = chrom_lengths[gene_info[gene].chrom]

        drugs = [d for d in gene_dict[gene] if d!=""]
        if len(drugs)==0:
            drugs = "None"
        else:
            drugs = ",".join(sorted(list(drugs)))
        lines.append([
            gene_info[gene].chrom,
            str(genome_start),
            str(genome_end),
            gene_info[gene].gene_id,
            gene_info[gene].name,
            drugs
        ])
    with open(outfile,"w") as O:
        for line in sorted(lines,key=lambda x: (x[0],int(x[1]))):
            O.write("%s\n" %"\t".join(line))

def assign_gene_to_amplicon(genes,chrom,start,end):
    l = []
    for g in genes.values():
        if g.chrom!=chrom: continue
        overlap = set(range(g.feature_start,g.feature_end)).intersection(set(range(int(start),int(end))))
        if overlap:
            l.append((g.gene_id,g.name,len(overlap)))
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
    ref = Fasta(ref_seq)
    with open(outfile,"w") as O:
        for chrom,start,end,amplicon_name in ref.get_amplicons(primer_file):
            locus_tag,gene_name = assign_gene_to_amplicon(genes,chrom,start,end)
            drugs = ",".join(assign_amplicon_drugs(db,chrom,start,end))
            if drugs=="":
                drugs = "None"
            O.write(f"{chrom}\t{start}\t{end}\t{locus_tag}\t{gene_name}\t{drugs}\t{amplicon_name}\n")

def so_term_in_mutation(mutation: str) -> bool:
    for term in supported_so_terms:
        if term in mutation:
            return True
    return False

def get_snpeff_formated_mutation_list(hgvs_variants,ref,gff,snpEffDB,db_dir):
    logging.debug("Converting HGVS to snpEff format")
    genes = load_gff(gff)
    refseq = FastaFile(ref)
    converted_mutations = {}
    so_term_rows = [r for r in hgvs_variants if so_term_in_mutation(r['Mutation'])]
    hgvs_variants = [r for r in hgvs_variants if not so_term_in_mutation(r['Mutation'])]

    converted_mutations = verify_mutation_list(hgvs_variants,genes,refseq, snpEffDB,db_dir)
    for row in so_term_rows:
        converted_mutations[(row["Gene"],row['Mutation'])] = (row['Gene'],row['Mutation'])

    return converted_mutations
    
def get_exon_to_aa_coords(exons):
    converter = []
    offset = 0
    if exons[0].strand == "+":
        for e in exons:
            aa = [math.floor((x - e.phase)/3) + 1 + offset for x in range(e.start - e.start  , e.end - e.start +1)]
            genome = range(e.start,e.end+1)
            converter.extend(list(zip(genome,aa)))
            offset = aa[-1]
    else:
        for e in exons[::-1]:
            aa = [math.floor((x - e.phase)/3) + 1 + offset for x in range(e.start - e.start  , e.end - e.start +1)]
            genome = range(e.end,e.start-1,-1)
            converter.extend(list(zip(genome,aa)))
            offset = aa[-1]
    return converter

def get_aa2genome_coords(exons):
    converter = get_exon_to_aa_coords(exons)
    aa2genome = defaultdict(list)
    for g,a in converter:
        aa2genome[a].append(g)
    return aa2genome

def get_genome_position(gene_object,change):
    g = gene_object
    for term in supported_so_terms:
        if term in change:
            return None
    if "any_missense_codon" in change:
        codon = int(change.replace("any_missense_codon_",""))
        change = f"p.Xyz{codon}Xyz"
 
    if change[0]=="p":
        aa2genome = get_aa2genome_coords(g.transcripts[0].exons)


    r = re.search("p.[A-Za-z]+([0-9]+)",change)
    if r:
        codon = int(r.group(1))
        return aa2genome[codon]
    
    # p.1?
    r = re.search(r'p.1\?',change)
    if r:
        codon = 1
        return aa2genome[codon]
    
    r = re.search("c.(-[0-9]+)[ACGT]+>[ACGT]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos
            return [p]
        else:
            p = g.start - pos
            return [p]
    r = re.search("n.(-?[0-9]+)[ACGT]+>[ACGT]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos -1
            return [p]

    # c.-30_-29insGCG
    r = re.search(r'[nc].([\-\*0-9]+)_([\-\*0-9]+)ins[A-Z]+',change)
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

    # n.211_212delGC
    r = re.search(r'[nc].([\-\*0-9]+)_([\-\*0-9]+)del[A-Z]*',change)
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

    # c.-37delT
    r = re.search(r'[nc].([\-\*0-9]+)del[A-Z]+',change)
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

    # n.1089dupC
    r = re.search(r'[nc].([\-0-9]+)dup[A-Z]+',change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]
    
    # c.134_135dupTC
    r = re.search(r'[nc].([\-0-9]+)_([\-0-9]+)dup[A-Z]+',change)
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

    

    r = re.search("[c].([0-9]+)([ACGT])>([ACGT])",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]
    
    r = re.search("g.([0-9]+)",change)
    if r:
        pos = int(r.group(1))
        # todo - check if this is correct add in chromosome
        return [pos]

    quit(f"Don't know how to handle {str(vars(g))} {change}")

def get_chrom_sizes(ref: FastaFile) -> dict:
    return dict(zip(ref.references, ref.lengths))

def match_ref_chrom_names(source: str,target: str) -> dict:
    logging.debug("Matching chromosome names")
    source_fa = FastaFile(source)
    source_fa_size = get_chrom_sizes(source_fa)
    target_fa = FastaFile(target)
    target_fa_size = get_chrom_sizes(target_fa)
    conversion = {}
    for s in target_fa.references:
        tlen = target_fa_size[s]
        tmp = [x[0] for x in source_fa_size.items() if x[1]==tlen]
        if len(tmp)==1:
            conversion[s] = tmp[0]
    return conversion

def replace_file_column(oldfilename,newfilename,column,conversion,sep="\t"):
    
    with open(newfilename,"w") as f:
        for line in open(oldfilename):
            row = line.strip().split(sep)
            if len(row)<column: continue
            for key,val in conversion.items():
                if row[column-1]==key:
                    row[column-1] = val
            f.write(sep.join(row)+"\n")


def create_db(args,extra_files = None):
    # set up snpeff folders if they don't exist
    create_snpeff_directories(args.db_dir)

    variables = json.load(open("variables.json"))    
    genome_file = "%s.fasta" % args.prefix
    gff_file = "%s.gff" % args.prefix
    bed_file = "%s.bed" % args.prefix
    json_file = "%s.dr.json" % args.prefix

    if os.path.isfile("snpEffectPredictor.bin"):
            snpeff_db_name = json.load(open("variables.json"))["snpEff_db"]
            load_snpEff_db("snpEffectPredictor.bin",snpeff_db_name,args.db_dir)
            
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
                row = l.strip().split("\t")
                if row[0] in chrom_conversion:
                    row[0] = chrom_conversion[row[0]]
                    O.write("\t".join(row)+"\n")        

    genes = load_gff(gff_file)
    gene_name2gene_id = {g.name:g.gene_id for g in genes}
    gene_name2gene_id.update({g.gene_id:g.gene_id for g in genes})
    gene_dict = {g.gene_id:g for g in genes}
    db = {}
    locus_tag_to_ann_dict = defaultdict(set)
    with open(args.prefix+".conversion.log","w") as L:
        if args.csv:
            hgvs_variants = [r for r in csv.DictReader(open(args.csv))]
            mutation_lookup = get_snpeff_formated_mutation_list(hgvs_variants,"genome.fasta","genome.gff",json.load(open("variables.json"))["snpEff_db"],args.db_dir)
            for row in csv.DictReader(open(args.csv)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                # annotation_info = {key:val for key,val in row.items() if key not in ["Gene","Mutation"]}
                # data looks like this: type=drug_resistance;drug=macrolides;literature=10.1038/s41467-021-25484-9
                annotation_info = {k:v for k,v in [x.split("=") for x in row["Info"].split(";")]}
                mut = mutation_lookup[(row["Gene"],row["Mutation"])][1]
                if mut!=row["Mutation"]:
                    L.write(f"Converted {row['Gene']} {row['Mutation']} to {mut}\n")
                if "drug" in annotation_info:
                    locus_tag_to_ann_dict[locus_tag].add(annotation_info["drug"])
                else:    
                    locus_tag_to_ann_dict[locus_tag].add(annotation_info["type"])
                if locus_tag not in db:
                    db[locus_tag] = {}
                if mut not in db[locus_tag]:
                    db[locus_tag][mut] = {"annotations":[]}

                tmp_annotation = annotation_info

                db[locus_tag][mut]["annotations"].append(tmp_annotation)
                db[locus_tag][mut]["genome_positions"] = get_genome_position(gene_dict[locus_tag],mut) if mut not in supported_so_terms else None
                db[locus_tag][mut]["chromosome"] = gene_dict[locus_tag].chrom
        

        gene_coordinates = {}
        if args.watchlist:
            for row in csv.DictReader(open(args.watchlist)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                locus_tag_to_ann_dict[locus_tag].add("")
                if row['Info']=="": continue
                info = {k:v for k,v in [x.split("=") for x in row["Info"].split(";")]}
                if 'drug' in info:
                    locus_tag_to_ann_dict[locus_tag].add(info['drug'])
                if 'start' in info and 'end' in info:
                    if locus_tag not in gene_coordinates:
                        gene_coordinates[locus_tag] = {
                            'chromosome':gene_dict[locus_tag].chrom,
                            'start':int(info['start']),
                            'end':int(info['end'])
                        }



        version = {"name":args.prefix}
        version.update(get_git_repo_info())

        version_obj = version.copy()
        for k,v in variables.items():
            if 'version' in k:
                version_obj[k] = v
        variables['version'] = version_obj
        # if 'db-schema-version' in variables:
        #     variables['version']['db-schema-version'] = variables['db-schema-version']

        json.dump(db,open(json_file,"w"),indent=4)
        
        
        for file in extra_files.values():
            if  isinstance(file,str):
                target = f"{args.prefix}.{file}"
                shutil.copyfile(file,target)
            else:
                target = f"{args.prefix}.{file['name']}"
                replace_file_column(file['name'],target,column=file['convert'],conversion=chrom_conversion)

        
        if "barcode" in extra_files:
            barcode_file = f"{args.prefix}.{extra_files['barcode']}"

            with open(barcode_file,"w") as O:
                for l in open("barcode.bed"):
                    if l[0]=="#": continue
                    row = l.strip().split("\t")
                    row[0] = chrom_conversion[row[0]]
                    O.write("\t".join(row)+"\n")
        
        if "amplicon_primers" in vars(args) and args.amplicon_primers:
            write_amplicon_bed(genome_file,gene_dict,db,args.amplicon_primers,bed_file)
            variables['amplicon'] = True
        else:
            write_bed(
                db=db,
                gene_dict=locus_tag_to_ann_dict,
                gene_info=gene_dict,
                ref_file=genome_file,
                outfile=bed_file,
                gene_coordinates=gene_coordinates if len(gene_coordinates)>0 else None,
            )
            variables['amplicon'] = False
        
                
        if list(chrom_conversion.keys())!=list(chrom_conversion.values()):
            variables["chromosome_conversion"] = {"target":list(chrom_conversion.keys()),"source":list(chrom_conversion.values())}
        variables_file = args.prefix+".variables.json"
        variables["files"] = {
            "ref": genome_file,
            "gff": gff_file,
            "bed": bed_file,
            "json_db": json_file,
            "variables": variables_file
        }
        if extra_files:
            for key,val in extra_files.items():
                if isinstance(val,str):
                    variables["files"][key] = f"{args.prefix}.{val}"
                else:
                    variables["files"][key] = f"{args.prefix}.{val['name']}"
                    
        json.dump(variables,open(variables_file,"w"))
        
        
        
        if args.load:
            load_db(variables_file,args.db_dir)

def index_ref(target):
    # pp.run_cmd(f"bwa index {target}")
    pp.run_cmd(f"samtools faidx {target}")
    tmp = target.replace(".fasta","")
    pp.run_cmd(f"samtools dict {target} -o {tmp}.dict")

def load_db(variables_file,db_dir,source_dir="."):
    variables = json.load(open(variables_file))

    if not os.path.isdir(db_dir):   
        os.mkdir(db_dir)

    for key,val in variables['files'].items():
        source = f"{source_dir}/{val}"
        target = f"{db_dir}/{val}"
        # logging.info(f"Copying file: {source} ---> {target}")
        shutil.copyfile(source,target)
        if key=="ref":
            index_ref(target)
    
    logging.info("[green]Sucessfully imported library[/]",extra={"markup":True})

def get_variable_file_name(db_dir:str,library_name:str):
    library_prefix = f"{db_dir}/{library_name}"
    return f"{library_prefix}.variables.json"

def check_db_files(variables):
    for key,val in variables.items():
        if key=="ref":
            if not os.path.isfile(val+".fai"):
                index_ref(val)

def is_db_path(string):
    if "/" in string:
        return True
    elif os.path.isfile(string+".variables.json"):
        return True
    return False

def check_db_exists(db_dir:str,db_name:str):
    db = get_db(db_dir=db_dir,db_name=db_name,verbose=False)
    if db is None:
        logging.error(f"DB {db_name} does not exist in the current directory or in {db_dir}")
        raise FileExistsError


def get_db(db_dir:str,db_name:str,verbose:bool=True):
    create_snpeff_directories(db_dir)
    if is_db_path(db_name):
        if "/" in db_name:
            db_dir = "/".join(db_name.split("/")[:-1])
            db_name = db_name.split("/")[-1]
            variable_file_name = f"{db_dir}/{db_name}.variables.json"
        else:
            db_dir = '.'
            variable_file_name = f"{db_dir}/{db_name}.variables.json"
    else:
        variable_file_name = os.path.join(db_dir,f"{db_name}.variables.json")
    
    if not os.path.isfile(variable_file_name):
        return None
    variables = json.load(open(variable_file_name))
    for key,val in variables['files'].items():
        if verbose:
            logging.info(f"Using {key} file: {db_dir}/{val}")
        if ".json" in val:
            variables[key] = json.load(open(f"{db_dir}/{val}"))
        elif key=='rules':
            variables[key] = [l.strip() for l in open(f'{db_dir}/{val}')]
        else:
            variables[key] = f"{db_dir}/{val}"
    
    check_db_files(variables)
    return variables    

def list_db(db_dir):
    if not os.path.isdir(db_dir):
        return []
    return [json.load(open(f"{db_dir}/{f}")) for f in os.listdir(db_dir) if f.endswith(".variables.json")]



def create_species_db(args: argparse.Namespace ,extra_files:dict = None, db_dir:str = None):
    variables = json.load(open("variables.json"))
    if not extra_files:
        extra_files = {}
    version = {
        "db-schema-version":variables['db-schema-version'],
        "name":args.prefix
    }
    for key,val in variables.items():
        if 'version' in key:
            version[key] = val
    
    version.update(get_git_repo_info())

    for file in extra_files.values():
        target = f"{args.prefix}.{file}"
        shutil.copyfile(file,target)

    variables_file = args.prefix+".variables.json"
    variables.update({
        "version": version,
        "files":{
            "variables": variables_file
        }
    })

    if extra_files:
        for key,val in extra_files.items():
            variables["files"][key] = f"{args.prefix}.{val}"
    json.dump(variables,open(variables_file,"w"))


    if not os.path.isdir(db_dir):   
        os.mkdir(db_dir)

    for key,val in variables['files'].items():
        target = f"{db_dir}/{val}"
        logging.debug(f"Copying file: {val} ---> {target}")
        shutil.copyfile(val,target)

def get_default_snpeff_dir():
    tmp = glob(f"{sys.base_prefix}/share/*snpeff*")
    if len(tmp)>0:
        return tmp[0]
    else: 
        return None

def get_default_snpeff_config():
    default_snpeff_dir = get_default_snpeff_dir()
    return f"{default_snpeff_dir}/snpEff.config"

def get_custom_snpeff_dir(db_dir):
    return f"{db_dir}/snpeff/"

def get_custom_snpeff_config(db_dir):
    return f"{get_custom_snpeff_dir(db_dir)}/snpEff.config"

def create_snpeff_directories(db_dir):
    default_snpeff_config = get_default_snpeff_config()
    custom_snpeff_dir = get_custom_snpeff_dir(db_dir)
    custom_snpeff_config = get_custom_snpeff_config(db_dir)
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)
    if not os.path.isdir(custom_snpeff_dir):
        os.mkdir(custom_snpeff_dir)
        os.mkdir(f"{custom_snpeff_dir}/data")
        with open(default_snpeff_config,"r") as INPUT, open(custom_snpeff_config,"w") as OUTPUT:
            for l in INPUT:
                OUTPUT.write(l)


def load_snpEff_db(bin_file: str,genome_name: str,db_dir:str):
    custom_snpeff_dir = get_custom_snpeff_dir(db_dir)
    custom_snpeff_config = get_custom_snpeff_config(db_dir)
    with open(custom_snpeff_config,"a") as F:
        F.write(f"{genome_name}.genome : {genome_name}\n")
    
    

    genome_dir = f"{custom_snpeff_dir}/data/{genome_name}"
    if not os.path.isdir(genome_dir):
        os.mkdir(genome_dir)
    shutil.copyfile(bin_file,f"{genome_dir}/{bin_file}")
