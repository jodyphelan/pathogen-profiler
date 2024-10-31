from collections import OrderedDict
from .utils import run_cmd, cmd_out
from .kmer import KmerDump
import os
import platform 
import logging
from .sourmash import SourmashSig
from .models import TargetQC, FastaQC
from typing import List
from .utils import shared_dict


class Fasta:
    """
    Class to represent fasta seuqnces in a python dict.

    Args:
        filename(str): Location of the fasta file

    Returns:
        fasta: A fasta class object
    """
    def __init__(self,filename):
        fa_dict = OrderedDict()
        seq_name = ""
        self.fa_file = filename
        for l in open(filename):
            line = l.rstrip()
            if line=="": continue
            if line.startswith(">"):
                seq_name = line[1:].split()[0]
                fa_dict[seq_name] = []
            else:
                fa_dict[seq_name].append(line)
        result = {}
        counter = 0
        sum_length = {}
        for seq in fa_dict:
            result[seq] = "".join(fa_dict[seq])
            result[seq] = result[seq].upper()
            sum_length[(counter+1,counter+len(result[seq]))] = seq
            counter = counter+len(result[seq])
        self.sum_length = sum_length
        self.fa_dict = result
    
    def align_to_ref(self,refseq,file_prefix):
        self.ref_aln = f"{file_prefix}.paf"
        run_cmd(f"minimap2 {refseq} {self.fa_file} --cs | sort -k6,6 -k8,8n > {self.ref_aln}")
        return self.ref_aln

    def get_amplicons(self,primer_file):
        bed = []
        for l in cmd_out(f"seqkit amplicon {self.fa_file} -p {primer_file} --bed"):
            row = l.strip().split()
            bed.append(tuple(row[:4]))
            
        return bed
    def get_kmer_counts(self,prefix,klen = 31,threads=1,max_mem=2, counter = "kmc"):
        shared_dict['kmer_counting'] = counter
        if counter=="kmc":
            if threads>32:
                threads = 32
            tmp_prefix = f"{prefix}_kmers"
            os.mkdir(tmp_prefix)
            bins = "-n128" if platform.system()=="Darwin" else ""
            run_cmd(f"kmc {bins} -m{max_mem} -t{threads} -k{klen} -ci1 -fm  {self.fa_file} {tmp_prefix} {tmp_prefix}")
            run_cmd(f"kmc_dump -ci1 {tmp_prefix} {prefix}.kmers.txt")
            run_cmd(f"rm -r {tmp_prefix}*")

            return KmerDump(f"{prefix}.kmers.txt",counter)
        elif counter=="dsk":
            max_mem = max_mem * 1000
            tmp_prefix = f"{prefix}_kmers"
            os.mkdir(tmp_prefix)
            run_cmd(f"dsk -file {self.fa_file} -abundance-min 1 -nb-cores {threads} -kmer-size {klen} -max-memory {max_mem} -out {tmp_prefix} -out-tmp {tmp_prefix}")
            run_cmd(f"dsk2ascii -file {tmp_prefix}.h5 -out {prefix}.kmers.txt")
            run_cmd(f"rm -r {tmp_prefix}*")

            return KmerDump(f"{prefix}.kmers.txt",counter)
    
    def sourmash_sketch(self,prefix,scaled=1000):
        logging.info("Sketching fasta")
        run_cmd(f"sourmash sketch dna -p abund,scaled={scaled} --merge {prefix} -o {prefix}.sig {self.fa_file}")
        return SourmashSig(f"{prefix}.sig",tmp_prefix=prefix)
    
    def get_n50(self):
        lengths = sorted([len(self.fa_dict[x]) for x in self.fa_dict])
        total = sum(lengths)
        half = total/2
        count = 0
        for l in lengths:
            count += l
            if count>=half:
                return l
            
    def get_fasta_qc(self) -> FastaQC:
        return FastaQC(
            num_sequences = len(self.fa_dict),
            num_bases = sum([len(self.fa_dict[x]) for x in self.fa_dict]),
            n50 = self.get_n50(),
            target_qc=[]
        )
    

class Paf:
    def __init__(self,filename: str):
        self.filename = filename

    def get_target_qc(self,bed_file: str) -> List[TargetQC]:
        results = []
        for l in cmd_out(f"cut -f6,8,9 {self.filename} | bedtools coverage -a {bed_file} -b -"):
            row = l.strip().split()
            results.append(
                TargetQC(
                    target = row[4],
                    percent_depth_pass=float(row[9])*100,
                    median_depth=int(float(row[9])),
                )
            )
        return results

    def get_ref_variants(self,refseq: str,sample_name: str,file_prefix: str) -> str:
        """
        Generate a vcf file of variants against a reference sequence from a paf file
        
        Arguments
        ---------
        refseq : str
            Reference sequence
        sample_name : str
            Sample name
        file_prefix : str
            Prefix for output file
        
        Returns
        -------
        str
            Filename of vcf file
        """
        self.refseq = refseq
        self.sample_name = sample_name
        self.file_prefix = file_prefix
        run_cmd("cat %(filename)s | paftools.js call -l 100 -L 100 -f %(refseq)s -s %(sample_name)s - | add_dummy_AD.py | bcftools view -Oz -o %(file_prefix)s.vcf.gz" % vars(self))
        return "%s.vcf.gz" % self.file_prefix
