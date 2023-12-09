import sys
import os.path
import subprocess
from collections import defaultdict
import random
import math
import re
import json
import csv
import subprocess as sp
from uuid import uuid4
import logging 
import pysam
from typing import List, NewType
from joblib import Parallel, delayed
from tqdm import tqdm

tmp_prefix = str(uuid4())

def get_tmp_file(prefix=None):
    """Get a temporary file"""
    if prefix:
        return "%s.%s" % (prefix,uuid4())
    else:
        return "%s" % uuid4()


def sanitize_region(region: str) -> str:
    """Replace : and - with _"""
    return region.replace(":","_").replace("-","_")

Region = NewType("Region",str)

def genome_job(cmd: str,region: Region) -> sp.CompletedProcess:
    """Run a command on a region of the genome"""
    region_safe = sanitize_region(region)
    cmd = cmd.format(region=region,region_safe=region_safe)
    out = run_cmd(cmd)
    return out


def get_genome_chunks(fasta: str,nchunks: int) -> List[Region]:
    """Split genome into n chunks"""
    genome = pysam.FastaFile(fasta)
    total_length = sum(genome.lengths)
    chunk_length = int(total_length/nchunks) + 10
    chunks = []
    chunk_start = 0
    chunk_end = 0
    for chrom in genome.references:
        while chunk_end < genome.get_reference_length(chrom):
            chunk_end += chunk_length
            if chunk_end > genome.get_reference_length(chrom):
                chunk_end = genome.get_reference_length(chrom)
            chunks.append([chrom,chunk_start,chunk_end])
            chunk_start = chunk_end
    regions = [f"{r[0]}:{r[1]+1}-{r[2]}" for r in chunks]
    return regions

def load_bed_regions(bed_file: str) -> List[Region]:
    """Load regions from a bed file"""
    regions = []
    with open(bed_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom,start,end = line.strip().split("\t")[:3]
            regions.append(f"{chrom}:{start}-{end}")
    return regions

def run_cmd_parallel_on_genome(cmd: str,genome: str,threads: int = 2,desc: str=None,bed_file: str=None) -> List[sp.CompletedProcess]:
    """Run a command in parallel in chunks of the genome"""
    logging.debug("Running command in parallel: %s" % cmd)
    if bed_file:
        regions = load_bed_regions(bed_file)
    else:
        regions = get_genome_chunks(genome,nchunks=threads)
    
    parallel = Parallel(n_jobs=threads, return_as="generator")
    desc = desc if desc else "Running command in parallel..."
    results = [r for r in tqdm(parallel(delayed(genome_job)(cmd,r) for r in regions),total=len(regions),desc=desc)]
    return results


def var_qc_test(var,min_depth,min_af,strand_support):
    """Test if a variant passes QC"""
    fail = False
    if min_depth!=None and var['depth']<min_depth:
        fail = True
    if min_af!=None and var['freq']<min_af:
        fail = True
    if strand_support!=None and var['forward_reads']!=None and var['forward_reads']<strand_support:
        fail = True
    if strand_support!=None and var['reverse_reads']!=None and var['reverse_reads']<strand_support:
        fail = True
    return fail

def sv_var_qc_test(var,min_depth,min_af,sv_len):
    
    fail = False
    if min_depth!=None and var['depth']<min_depth:
        fail = True
    if min_af!=None and var['freq']<min_af:
        fail = True
    if sv_len!=None and var["sv_len"]>sv_len:
        fail = True
    return fail

def filter_variant(var,filter_params):
    qc = "pass"
    if var['sv']==True:
        if sv_var_qc_test(var,filter_params["sv_depth_hard"],filter_params["sv_af_hard"],filter_params["sv_len_hard"]):
            qc = "hard_fail"
        elif sv_var_qc_test(var,filter_params["sv_depth_soft"],filter_params["sv_af_soft"],filter_params["sv_len_soft"]):
            qc = "soft_fail"
    else:
        if var_qc_test(var,filter_params["depth_hard"],filter_params["af_hard"],filter_params["strand_hard"]):
            qc = "hard_fail"
        elif var_qc_test(var,filter_params["depth_soft"],filter_params["af_soft"],filter_params["strand_soft"]):
            qc = "soft_fail"
    return qc


def stringify(l):
    return [str(x) for x in list(l)]

def parse_csv(filename):
    """Parses a CSV file into a dictionary using the first column as the key."""
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        key = reader.fieldnames[0]
        return {row[key]: row for row in reader}

def return_fields(obj,args,i=0):
    largs = args.split(".")
    if i+1>len(largs):
        return obj
    if largs[i] not in obj:
        return None
    sub_obj = obj[largs[i]]
    if isinstance(sub_obj,dict):
        return return_fields(sub_obj,args,i+1)
    elif isinstance(sub_obj,list):
        return [return_fields(x,args,i+1) for x in sub_obj]
    else:
        return sub_obj
        
def variable2string(var,quote=False):
    q = '"' if quote else ""
    if isinstance(var,float):
        return "%.3f" % var
    elif isinstance(var,dict):
        return "%s%s%s" % (q,",".join(list(var)),q)
    elif isinstance(var,list):
        return "%s%s%s" % (q,",".join(var),q)
    else:
        return "%s%s%s" % (q,str(var),q)

def dict_list2text(l,columns = None, mappings = None,sep="\t"):
    if mappings:
        headings = list(mappings)
    elif columns:
        headings = columns
    else:
        headings = list(l[0].keys())
    rows = []
    header = sep.join([mappings[x].title() if (mappings!=None and x in mappings) else x.title() for x in headings])
    for row in l:
        r = sep.join([variable2string(return_fields(row,x)) for x in headings])
        rows.append(r)
    str_rows = "\n".join(rows)
    out ="%s\n%s" % (header,str_rows)
    return out.strip()


def get_lt2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[3]] = None if row[5]=="None" else row[5].split(",") 
    return lt2drugs


def process_variants(results: dict,conf: dict,annotations: List[str]):
    variant_containers = {d:[] for d in annotations}
    variant_containers['other'] = []
    variant_containers['qc_fail'] = []
    for var in results['variants']:
        annotation_containers = {d:[a for a in var['annotation'] if a['type']==d] for d in annotations}
        qc  = filter_variant(var,conf["variant_filters"])
        if qc=="hard_fail":
            continue
        elif qc=="soft_fail":
            variant_containers['qc_fail'].append(var)
        else:
            assigned = False
            for a in annotations:
                if annotation_containers[a]:
                    assigned = True
                    variant_containers[a].append(var)
            if not assigned:
                variant_containers['other'].append(var)
            
    for a in annotations:
        results[a+"_variants"] = variant_containers[a]
    results['other_variants'] = variant_containers['other']
    results['qc_fail_variants'] = variant_containers['qc_fail']
    del results['variants']

    return results



def reformat_annotations(results,conf):
    #Chromosome      4998    Rv0005  -242
    lt2drugs = get_lt2drugs(conf["bed"])
    results["dr_variants"] = []
    results["other_variants"] = []
    results["qc_fail_variants"] = []
    for var in results["variants"]:
        drugs = tuple([x["drug"] for x in var.get("annotation",[]) if x["type"]=="drug_resistance"])
        if len(drugs)>0:
            tmp = var.copy()
            dr_ann = []
            other_ann = []
            while len(tmp["annotation"])>0:
                x = tmp["annotation"].pop()
                if x["type"]=="drug_resistance":
                    dr_ann.append(x)
                else:
                    other_ann.append(x)
            tmp["drugs"] = dr_ann
            tmp["annotation"] = other_ann
            tmp["gene_associated_drugs"] = lt2drugs[var["locus_tag"]]
            qc  = filter_variant(var,conf["variant_filters"])
            if qc=="hard_fail":
                continue
            elif qc=="soft_fail":
                results["qc_fail_variants"].append(tmp)
            else:
                results["dr_variants"].append(tmp)
        else:
            var["gene_associated_drugs"] = lt2drugs[var["locus_tag"]]
            qc  = filter_variant(var,conf["variant_filters"])
            if qc=="hard_fail":
                continue
            elif qc=="soft_fail":
                results["qc_fail_variants"].append(var)
            else:
                results["other_variants"].append(var)
    del results["variants"]
    return results

def get_genome_positions_from_db(db):
    genome_positions = defaultdict(set)
    for gene in db:
        for var in db[gene]:
            drugs = tuple([x["drug"] for x in db[gene][var]["annotations"] if x["type"]=="drug"])
            if len(drugs)==0:
                continue
            if db[gene][var]["genome_positions"]:
                for pos in db[gene][var]["genome_positions"]:
                    genome_positions[pos].add((gene,var,drugs))

    return genome_positions

def lt2genes(bed_file):
    #Chromosome      759310  763325  Rv0667  rpoB    rifampicin
    lt2gene = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2gene[row[3]] = row[4]
    return lt2gene


def reformat_missing_genome_pos(positions,conf):
    lt2gene = lt2genes(conf["bed"])
    dr_associated_genome_pos = get_genome_positions_from_db(conf["json_db"])
    new_results = []
    for pos in positions:
        if pos in dr_associated_genome_pos:
            tmp = dr_associated_genome_pos[pos]
            gene = list(tmp)[0][0]
            variants = ",".join([x[1] for x in tmp])
            drugs = ",".join(set(unlist([x[2] for x in tmp])))
            new_results.append({"position":pos,"locus_tag":gene, "gene": lt2gene[gene], "variants": variants, "drugs":drugs})
    return new_results

def select_most_relevant_csq(csqs):
    rank = ["transcript_ablation","exon_loss_variant","frameshift_variant","large_deletion","start_lost","disruptive_inframe_deletion","disruptive_inframe_insertion","stop_gained","stop_lost","conservative_inframe_deletion","conservative_inframe_insertion","initiator_codon_variant","missense_variant","non_coding_transcript_exon_variant","upstream_gene_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_variant","stop_retained_variant","splice_region_variant","synonymous_variant"]

    ranked_csq = sorted(csqs,key=lambda x: min([rank.index(y) if y in rank else 999 for y in x['type'].split("&")]))
    return ranked_csq[0]

def set_change(var):
    protein_csqs = ["missense_variant","stop_gained"]
    var["change"] = var["protein_change"] if var["type"] in protein_csqs else var["nucleotide_change"]
    return var

def annotation_has_drug_type(ann):
    return any([x["type"]=="drug" for x in ann])

def select_csq(dict_list):
    for d in dict_list:
        annotated_csq = []
        for csq in d["consequences"]:
            if "annotation" in csq:
                annotated_csq.append(csq)
        if len(annotated_csq)==0:
            csq = select_most_relevant_csq(d["consequences"])
            alternate_consequences = [json.dumps(x) for x in d["consequences"]]
            alternate_consequences.remove(json.dumps(csq))
            alternate_consequences = [json.loads(x) for x in alternate_consequences]
        elif len(annotated_csq)==1:
            csq = annotated_csq[0]
            alternate_consequences = []
        else:
            chosen_annotation = None
            for csq in annotated_csq:
                if annotation_has_drug_type(csq["annotation"]):
                    chosen_annotation = csq
                    break
            if chosen_annotation:
                csq = chosen_annotation
            else:
                csq = annotated_csq[0]
            alternate_consequences = [json.dumps(x) for x in d["consequences"]]
            alternate_consequences.remove(json.dumps(csq))
            alternate_consequences = [json.loads(x) for x in alternate_consequences]
        del d["consequences"]
        d.update(csq)
        d["alternate_consequences"] = alternate_consequences
        d = set_change(d)
    return dict_list

def dict_list_add_genes(dict_list,conf,key="gene_id"):
    lt2gene = {}
    for l in open(conf["bed"]):
        row = l.rstrip().split()
        lt2gene[row[3]] = row[4]
    for d in dict_list:
        d["locus_tag"] = d[key]
        d["gene"] = lt2gene[d[key]]
        del d[key]
        if "gene_name" in d:
            del d["gene_name"]
    return dict_list



def iupac(n):
    tmp = {
        "A":["A"],
        "C":["C"],
        "G":["G"],
        "T":["T"],
        "R":["A","G"],
        "Y":["C","T"],
        "S":["G","C"],
        "W":["A","T"],
        "K":["G","T"],
        "M":["A","C"],
        "B":["C","G","T"],
        "D":["A","G","T"],
        "H":["A","C","T"],
        "V":["A","C","G"],
        "N":["A","C","G","T"]
    }
    return tmp[n]

def unlist(t):
    return [item for sublist in t for item in sublist]
    

def get_seqs_from_bam(bamfile):
    seqs = []
    for l in cmd_out("samtools view %s -H | grep ^@SQ" % bamfile):
        row = l.rstrip().split()
        seqs.append(row[1].replace("SN:",""))
    return seqs


def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])


def stdev(arr):
    mean = sum(arr)/len(arr)
    return math.sqrt(sum([(x-mean)**2 for x in arr])/len(arr))


def add_arguments_to_self(self,args: dict) -> None:
    # Function to add arguments to class instance
    for x in args:
        if x == "self":
            continue
        vars(self)[x] = args[x]
    if "kwargs" in args:
        for x in args["kwargs"]:
            vars(self)[x] = args["kwargs"][x]



def run_cmd(cmd: str, desc=None, log: str=None) -> sp.CompletedProcess:
    if desc:
        logging.info(desc)
    programs = set([x.strip().split()[0] for x in re.split("[|&;]",cmd) if x!=""])
    missing = [p for p in programs if which(p)==False]
    if len(missing)>0:
        raise ValueError("Cant find programs: %s\n" % (", ".join(missing)))
    logging.debug(f"Running command: {cmd}")
    cmd = "/bin/bash -c set -o pipefail; " + cmd
    output = open(log,"w") if log else sp.PIPE
    result = sp.run(cmd,shell=True,stderr=output,stdout=output)
    if result.returncode != 0:
        logging.error(result.stderr.decode("utf-8"))
        raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,result.stderr.decode()))
    return result

def cmd_out(cmd: str) -> str:
    filename = str(uuid4())
    cmd = f"{cmd} > {filename}"
    run_cmd(cmd)
    for line in open(filename):
        yield line.strip()
    os.remove(filename)

def _cmd_out(cmd,verbose=1):
    cmd = "set -u pipefail; " + cmd
    if verbose==2:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/stderr","w")
    elif verbose==1:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/null","w")
    else:
        stderr = open("/dev/null","w")
    try:
        res = subprocess.Popen(cmd,shell=True,stderr = stderr,stdout=subprocess.PIPE)
        for l in res.stdout:
            yield l.decode().rstrip()
    except:
        logging.error("Command Failed! Please Check!")
        raise Exception
    stderr.close()

def log(msg,ext=False):
    sys.stderr.write("\n"+str(msg)+"\n")
    if ext:
        exit(1)


def load_bed(filename,columns,key1,key2=None,intasint=False):
    results = defaultdict(lambda: defaultdict(tuple))
    for l in open(filename):
        row = l.rstrip().split("\t")
        if l[0]=="#":
            header = row
            continue
        if key2:
            if max(columns+[key1,key2])>len(row):
                logging.error("Can't access a column in BED file. The largest column specified is too big",True)
            if key2==2 or key2==3:
                results[row[key1-1]][int(row[key2-1])] = tuple([row[int(x)-1] for x in columns])
            else:
                results[row[key1-1]][row[key2-1]] = tuple([row[int(x)-1] for x in columns])
        else:
            if max(columns+[key1])>len(row):
                logging.error("Can't access a column in BED file. The largest column specified is too big",True)
            results[row[key1-1]]= tuple([row[int(x)-1] for x in columns])
    return results


def filecheck(filename):
    """
    Check if file is there and quit if it isn't
    """
    if filename=="/dev/null":
        return filename
    elif not os.path.isfile(filename):
        logging.error("Can't find %s\n" % filename)
        exit(1)
    else:
        return filename


def nofile(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isfile(filename):
        return True
    else:
        return False

def nofolder(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isdir(filename):
        return True
    else:
        return False

def create_seq_dict(ref):
    ref_prefix = ref.replace(".fasta","").replace(".fa","")
    if nofile(f"{ref_prefix}.dict"):
        run_cmd("samtools dict %s -o {ref_prefix}.dict" % ref)

def bowtie_index(ref):
    if nofile("%s.1.bt2"%ref):
        cmd = "bowtie2-build %s %s" % (ref,ref)
        run_cmd(cmd)

def bwa2_index(ref):
    """
    Create BWA index for a reference
    """
    if nofile("%s.bwt.2bit.64"%ref):
        cmd = "bwa-mem2 index %s" % ref
        run_cmd(cmd)

def bwa_index(ref):
    """
    Create BWA index for a reference
    """
    if nofile("%s.bwt"%ref):
        cmd = "bwa index %s" % ref
        run_cmd(cmd)

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True
    return False


# def run_cmd(cmd,verbose=1,target=None,terminate_on_error=True):
#     """
#     Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
#     """
#     programs = set([x.strip().split()[0] for x in re.split("[|&;]",cmd) if x!=""])
#     missing = [p for p in programs if which(p)==False]
#     if len(missing)>0:
#         raise ValueError("Cant find programs: %s\n" % (", ".join(missing)))
#     if target and filecheck(target): return True
#     cmd = "set -u pipefail; " + cmd
#     if verbose>0:
#         sys.stderr.write("\nRunning command:\n%s\n" % cmd)

#     p = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     stdout,stderr = p.communicate()

#     if terminate_on_error is True and p.returncode!=0:
#         raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,stderr.decode()))

#     if verbose>1:
#         sys.stdout.write(stdout)
#         sys.stderr.write(stderr)

#     return (stdout.decode(),stderr.decode())

def index_bam(bamfile,threads=1,overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "samtools index -@ %s %s" % (threads,bamfile)
    bam_or_cram = "cram" if bamfile[-4:]=="cram" else "bam"
    suffix = ".bai" if bam_or_cram=="bam" else ".crai"
    if filecheck(bamfile):
        if nofile(bamfile+suffix):
            run_cmd(cmd)
        elif os.path.getmtime(bamfile+suffix)<os.path.getmtime(bamfile) or overwrite:
            run_cmd(cmd)

def index_bcf(bcffile,threads=1,overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "bcftools index --threads %s -f %s" % (threads,bcffile)
    if filecheck(bcffile):
        if nofile(bcffile+".csi"):
            run_cmd(cmd)
        elif os.path.getmtime(bcffile+".csi")<os.path.getmtime(bcffile) or overwrite:
            run_cmd(cmd)

def tabix(bcffile,threads=1,overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "bcftools index --threads %s -ft %s" % (threads,bcffile)
    if filecheck(bcffile):
        if nofile(bcffile+".tbi"):
            run_cmd(cmd)
        elif os.path.getmtime(bcffile+".tbi")<os.path.getmtime(bcffile) or overwrite:
            run_cmd(cmd)

def rm_files(x):
    """
    Remove a files in a list format
    """
    for f in x:
        if os.path.isfile(f):
            logging.debug("Removing %s" % f)
            os.remove(f)


