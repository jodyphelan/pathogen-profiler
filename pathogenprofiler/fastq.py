from __future__ import division
from .bam import Bam
from .utils import filecheck, add_arguments_to_self, run_cmd,bwa_index,bwa_meme_index,bwa2_index,bowtie_index,cmd_out, shared_dict
import os
from collections import defaultdict
from .kmer import KmerDump
import platform
import logging
from .sourmash import SourmashSig
from .models import FastqQC



class Fastq:
    """
    Class to hold fastq file and methods.
    Methods include trimming and mapping to a reference genome
    """
    def __init__(self,r1,r2=None,r3=None):
        """
        r1 = Forward reads
        r2 = Reverse reads (optional)
        r3 = Unpaired reads (optional)
        """
        add_arguments_to_self(self, locals())
        # Work out if it is paired end sequencing
        self.paired = True if (r1 and r2) else False
        filecheck(r1)
        self.files = [r1]
        if self.paired:
            filecheck(r2)
            self.files.append(r2)
        if r3:
            filecheck(r3)
            self.files.append(r3)

    def trim(self, prefix, threads=1):
        """Perform trimming"""
        logging.info("Trimming reads")
        shared_dict['trimming'] = 'trimmomatic'
        add_arguments_to_self(self, locals())
        if self.paired:
            run_cmd("trimmomatic PE -threads %(threads)s -phred33 %(r1)s %(r2)s -baseout %(prefix)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(self))
            run_cmd("cat %(prefix)s_1U %(prefix)s_2U > %(prefix)s_TU" % vars(self))
            run_cmd("rm %(prefix)s_1U %(prefix)s_2U" % vars(self))
            return Fastq("%(prefix)s_1P" % vars(self), "%(prefix)s_2P" % vars(self), "%(prefix)s_TU" % vars(self))
        else:
            run_cmd("trimmomatic SE -threads %(threads)s -phred33 %(r1)s %(prefix)s_TU LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(self))
            return Fastq("%(prefix)s_TU" % vars(self))

    def map_to_ref(self, ref_file, prefix, sample_name, aligner, platform, threads=1,markdup=True, max_mem="768M"):
        """Mapping to a reference genome"""
        logging.info("Mapping to reference genome")
        add_arguments_to_self(self, locals())
        self.aligner = aligner.lower()
        accepted_aligners = ["bwa","bwa-meme","bwa-mem2","bowtie2","minimap2"]
        if self.aligner not in accepted_aligners:
            quit("ERROR: %s not in accepted aligners\n" % aligner)

        self.platform = platform.lower()
        accepted_platforms = ["illumina","nanopore","pacbio"]
        if self.platform not in accepted_platforms:
            quit("ERROR: %s not in accepted platforms\n" % platform)

        if self.aligner=="minimap2":
            pass
        else:
            {"bwa":bwa_index,"bwa-meme":bwa_meme_index,"bwa-mem2":bwa2_index,"bowtie2":bowtie_index}[self.aligner](ref_file)

        self.bwa_prefix = "bwa mem -t %(threads)s -K 10000000 -c 100 -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -M -T 50" % vars(self)
        self.bwa_meme_prefix = "bwa-meme mem -t %(threads)s -K 10000000 -c 100 -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -M -T 50" % vars(self)
        self.bwa2_prefix = "bwa-mem2 mem -t %(threads)s -c 100 -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -M -T 50" % vars(self)
        self.bowtie2_prefix = "bowtie2 -p %(threads)s --rg-id '%(sample_name)s' --rg 'SM:%(sample_name)s' --rg 'PL:%(platform)s'" % vars(self)
        self.minimap2_prefix = "minimap2 --MD -t %(threads)s -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -a" % vars(self)
        self.bam_file = "%s.bam" % self.prefix
        self.bam_single_file = "%s.single.bam" % self.prefix
        self.bam_pair_file = "%s.pair.bam" % self.prefix
        self.bam_unsort_file = "%s.unsort.bam" % self.prefix
        if self.platform == "nanopore":
            run_cmd("%(minimap2_prefix)s -x map-ont  %(ref_file)s %(r1)s | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
        elif self.platform == "pacbio":
            run_cmd("%(minimap2_prefix)s -x map-pb  %(ref_file)s %(r1)s | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
        else:
            if aligner=="bwa" and self.paired:
                run_cmd("%(bwa_prefix)s %(ref_file)s %(r1)s %(r2)s | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                if self.r3:
                    run_cmd("%(bwa_prefix)s %(ref_file)s %(r3)s | samtools sort -@ %(threads)s -o %(bam_single_file)s -" % vars(self))
            elif aligner=="bwa" and not self.paired:
                run_cmd("%(bwa_prefix)s %(ref_file)s %(r1)s | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
            elif aligner=="bwa-meme" and self.paired:
                run_cmd("%(bwa_meme_prefix)s %(ref_file)s %(r1)s %(r2)s | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                if self.r3:
                    run_cmd("%(bwa_meme_prefix)s %(ref_file)s %(r3)s | samtools sort -@ %(threads)s -o %(bam_single_file)s -" % vars(self))
            elif aligner=="bwa-meme" and not self.paired:
                run_cmd("%(bwa_meme_prefix)s %(ref_file)s %(r1)s | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
            elif aligner=="bwa-mem2" and self.paired:
                run_cmd("%(bwa2_prefix)s %(ref_file)s %(r1)s %(r2)s | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                if self.r3:
                    run_cmd("%(bwa2_prefix)s %(ref_file)s %(r3)s | samtools sort -@ %(threads)s -o %(bam_single_file)s -" % vars(self))
            elif aligner=="bwa-mem2" and not self.paired:
                run_cmd("%(bwa2_prefix)s %(ref_file)s %(r1)s | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
            elif aligner=="bowtie2" and self.paired:
                run_cmd("%(bowtie2_prefix)s -x %(ref_file)s -1 %(r1)s -2 %(r2)s | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                if self.r3:
                    run_cmd("%(bowtie2_prefix)s -x %(ref_file)s -U %(r3)s | samtools sort -@ %(threads)s -o %(bam_single_file)s" % vars(self))
            elif aligner=="bowtie2" and not self.paired:
                run_cmd("%(bowtie2_prefix)s  -x %(ref_file)s -1 %(r1)s | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
            elif aligner=="minimap2" and self.paired:
                run_cmd("%(minimap2_prefix)s -ax sr %(ref_file)s %(r1)s %(r2)s | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                if self.r3:
                    run_cmd("%(minimap2_prefix)s -ax sr %(ref_file)s %(r3)s| samtools sort -@ %(threads)s -o %(bam_single_file)s -" % vars(self))
            elif aligner=="minimap2" and not self.paired:
                run_cmd("%(minimap2_prefix)s -ax sr %(ref_file)s %(r1)s| samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))

            if self.paired:
                if self.r3:
                    run_cmd("samtools merge -@ %(threads)s -f %(bam_unsort_file)s %(bam_pair_file)s %(bam_single_file)s" % vars(self))
                else:
                    self.bam_unsort_file = self.bam_pair_file
                # run_cmd("samtools sort -@ %(threads)s -o %(bam_file)s %(bam_unsort_file)s" % vars(self))
                if markdup:
                    run_cmd("samtools sort -m %(max_mem)s -n -@ %(threads)s  %(bam_unsort_file)s | samtools fixmate -@ %(threads)s -m - - | samtools sort -m %(max_mem)s -@ %(threads)s - | samtools markdup -@ %(threads)s - %(bam_file)s" % vars(self))
                else:
                    run_cmd("samtools sort -m %(max_mem)s -@ %(threads)s -o %(bam_file)s %(bam_unsort_file)s" % vars(self))
                if self.r3:
                    run_cmd("rm %(bam_single_file)s %(bam_pair_file)s %(bam_unsort_file)s" % vars(self))
                else:
                    run_cmd("rm %(bam_pair_file)s" % vars(self))
        shared_dict['mapping'] = aligner
        return Bam(self.bam_file,self.prefix,self.platform,threads=threads)
    
    def get_kmer_counts(self,prefix,klen = 31,threads=1,max_mem=8,counter = "kmc"):
        logging.info("Counting kmers")
        shared_dict['kmer_counting'] = counter
        if counter=="kmc":
            if threads>32:
                threads = 32
            tmp_prefix = f"{prefix}_kmers"
            tmp_file_list = f"{prefix}.kmc.list"
            os.mkdir(tmp_prefix)
            with open(tmp_file_list,"w") as O:
                O.write("\n".join(self.files))
            bins = "-n128" if platform.system()=="Darwin" else ""
            run_cmd(f"kmc {bins} -m{max_mem} -t{threads} -k{klen} @{tmp_file_list} {tmp_prefix} {tmp_prefix}")
            run_cmd(f"kmc_dump {tmp_prefix} {prefix}.kmers.txt")
            run_cmd(f"rm -r {tmp_prefix}*")

            return KmerDump(f"{prefix}.kmers.txt",counter)
        elif counter=="dsk":
            max_mem = max_mem * 1000
            tmp_prefix = f"{prefix}_kmers"
            os.mkdir(tmp_prefix)
            r2 = f"-file {self.r2}" if self.r2 else ""
            run_cmd(f"dsk -file {self.r1} {r2} -abundance-min 2 -nb-cores {threads} -kmer-size {klen} -max-memory {max_mem} -out {tmp_prefix} -out-tmp {tmp_prefix}")
            run_cmd(f"dsk2ascii -file {tmp_prefix}.h5 -out {prefix}.kmers.txt")
            run_cmd(f"rm -r {tmp_prefix}*")

            return KmerDump(f"{prefix}.kmers.txt",counter)
    
    def sourmash_sketch(self,prefix,scaled=1000):
        logging.info("Sketching reads")
        read1 = self.r1
        read2 = self.r2 if self.r2 else ""
        run_cmd(f"sourmash sketch dna -p abund,scaled={scaled} --merge {prefix} -o {prefix}.sig {read1} {read2}")
        return SourmashSig(f"{prefix}.sig",tmp_prefix=prefix)
    
    def get_qc(self):
        """Get quality control metrics"""
        header = None
        result = defaultdict(int)
        for l in cmd_out('seqkit stats -T %s' % "\t".join(self.files)):
            row = l.strip().split()
            if row[0]=='file': 
                header = row
                continue
            r = dict(zip(header,row))
            for column in ('num_seqs','sum_len'):
                result[column] += int(r[column])
        return FastqQC(
            num_bases=result['sum_len'],
            num_sequences=result['num_seqs']
        )