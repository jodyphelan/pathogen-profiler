from glob import glob
from .kmer import KmerDump
from .utils import TempFilePrefix, load_bed_regions, get_genome_chunks, add_arguments_to_self, run_cmd, cmd_out, filecheck, index_bam, run_cmd_parallel_on_genome, load_bed
from .vcf import Vcf
from collections import defaultdict
from uuid import uuid4
import os
import platform 
import statistics as stats
import logging 
from pysam import FastaFile
from .models import BamQC, TargetQC, GenomePositionDepth, GenomePosition
from typing import List, Optional
from .utils import shared_dict




class Bam:
    """
    A class to perform operations on BAM files such as SNP calling and QC
    """
    def __init__(
        self,
        bam_file: str,
        prefix: str,
        platform: str,
        threads: int = 1
    ):
        self.bam_file = bam_file
        self.prefix = prefix
        self.platform = platform
        self.threads = threads


        logging.debug("Creating Bam object with %s" % (bam_file))
        filecheck(self.bam_file)
        index_bam(bam_file,threads=threads)
        self.filetype = "cram" if bam_file[-5:]==".cram" else "bam"
        for l in cmd_out("samtools view -H %s" % (bam_file)):
            if l[:3]=="@RG":
                row = l.strip().split("\t")
                for r in row:
                    if r.startswith("SM:"):
                        self.bam_sample_name = r.replace("SM:","")
    def calculate_bed_depth(self,bed_file: str) -> List[GenomePositionDepth]:
        """
        Calculate depth of BAM file in regions specified by a BED file
        
        Arguments
        ---------
        bed_file : str
            BED file containing regions to calculate depth for
        
        Returns
        -------
        List of GenomePositionDepth objects
        """
        logging.info("Calculating depth in regions")
        self.bed_file = bed_file
        position_depth = defaultdict(list)
        for r, _ in load_bed(self.bed_file).items():
                for p in r.iter_positions():
                    position_depth[p] = 0


        for l in cmd_out("samtools view -Mb -L %(bed_file)s %(bam_file)s | samtools depth - " % vars(self)):
            row = l.strip().split()
            p = GenomePosition(chrom=row[0],pos=int(row[1]))
            position_depth[p] = int(row[2])

        self.position_depth = [GenomePositionDepth(chrom=p.chrom,pos=p.pos,depth=v) for p,v in position_depth.items()]
        return self.position_depth

    def run_delly(self,ref_file: str,bed_file: str) -> Vcf:
        """
        Method to run delly and extract variants overlapping with a BED file
        
        Arguments
        ---------
        bed_file : str
            BED file to extract variants from
        
        Returns
        -------
        Vcf object
        """
        logging.info("Running delly")
        self.bed_file = bed_file
        self.ref_file = ref_file
        if self.platform=="illumina":
            run_cmd("delly call -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self))            
        else:
            run_cmd("delly lr -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self))

        run_cmd("bcftools view -c 2 %(prefix)s.delly.bcf | bcftools view -e '(INFO/END-POS)>=100000' -Oz -o %(prefix)s.delly.vcf.gz" % vars(self))
        run_cmd("bcftools index %(prefix)s.delly.vcf.gz" % vars(self))
        run_cmd("bcftools view -R %(bed_file)s %(prefix)s.delly.vcf.gz -Oz -o %(prefix)s.delly.targets.vcf.gz" % vars(self))
        return Vcf("%(prefix)s.delly.targets.vcf.gz" % vars(self))

    def call_variants(
        self,
        ref_file: str,
        caller: str,
        filters: dict,
        bed_file: Optional[str] = None,
        threads: int = 1,
        calling_params: Optional[str] = None, 
        samclip: bool = False
    ) -> Vcf:
        from .variant_calling import VariantCaller
        subclasses = {cls.__software__:cls for cls in VariantCaller.__subclasses__()}
        chosen_class = subclasses[caller]
        caller = chosen_class(
            ref_file=ref_file,
            bam_file=self.bam_file,
            prefix=self.prefix,
            bed_file=bed_file,
            threads=threads,
            samclip=samclip,
            platform=self.platform,
            calling_params=calling_params,
            filters=filters
        )
        return caller.call_variants()
    
    
    
    def __call_variants__deprecated(
        self,
        ref_file: str,
        caller: str,
        filters: str,
        bed_file: Optional[str] = None,
        threads: int = 1,
        calling_params: Optional[str] = None, 
        samclip: bool = False
    ) -> Vcf:
        """Method to run variant calling"""
        add_arguments_to_self(self, locals())
        filecheck(ref_file)
        self.caller = caller.lower()
        self.af_hard = filters['af_hard'] if 'af_hard' in filters else 0
        # Set up final vcf file name
        # Make the windows for parallel calling based on chunking the whole
        # genome or by providing a bed file
        with TempFilePrefix() as tmp:
            self.temp_file_prefix = tmp
            if bed_file:
                self.vcf_file = "%s.short_variants.targets.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)
                genome_chunks = [r.safe for r in load_bed_regions(bed_file)]
            else:
                self.vcf_file = "%s.short_variants.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)
                genome_chunks = [r.safe for r in get_genome_chunks(self.ref_file,self.threads)]

            if self.platform in ("illumina"):
                self.samclip_cmd = "| samclip --ref %(ref_file)s" % vars(self) if samclip else ""
            else:
                self.samclip_cmd = ""
            
            # Run through different options. 

            # Nanopore
            if self.platform == "nanopore" and self.caller=="bcftools":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD,ADF,ADR -r {region} %(bam_file)s | bcftools call -mv | annotate_maaf.py | bcftools +fill-tags | bcftools norm -f %(ref_file)s | bcftools filter -e 'IMF < 0.7' -S 0 -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform == "nanopore" and self.caller=="freebayes":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "freebayes -f %(ref_file)s -F %(af_hard)s -r {region} --haplotype-length -1 %(calling_params)s %(bam_file)s | annotate_maaf.py  | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform=="nanopore" and self.caller == "pilon":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = """pilon --genome %(ref_file)s --targets {region} %(calling_params)s --nanopore %(bam_file)s --variant --output %(temp_file_prefix)s.{region_safe} && bcftools view -i 'AF>0' %(temp_file_prefix)s.{region_safe}.vcf | fix_pilon_headers.py --sample %(bam_sample_name)s| add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz""" % vars(self)
            
            # Pacbio
            elif self.platform == "pacbio" and self.caller=="freebayes":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "freebayes -f %(ref_file)s -r {region} --haplotype-length -1 %(calling_params)s %(bam_file)s | annotate_maaf.py | bcftools +fill-tags | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform=="pacbio" and self.caller == "pilon":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = """pilon --genome %(ref_file)s --targets {region} %(calling_params)s --pacbio %(bam_file)s --variant --output %(temp_file_prefix)s.{region_safe} && bcftools view -i 'AF>0' %(temp_file_prefix)s.{region_safe}.vcf | fix_pilon_headers.py --sample %(bam_sample_name)s| add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz""" % vars(self)
            
            # Illumina
            elif self.platform=="illumina" and self.caller == "bcftools":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD,ADF,ADR -r {region} %(temp_file_prefix)s.{region_safe}.tmp.bam | bcftools call -mv | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform=="illumina" and self.caller == "gatk":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && gatk HaplotypeCaller -R %(ref_file)s -I %(temp_file_prefix)s.{region_safe}.tmp.bam -O /dev/stdout -L {region} %(calling_params)s -OVI false | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform=="illumina" and self.caller == "freebayes":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && freebayes -f %(ref_file)s -r {region} --haplotype-length -1 %(calling_params)s %(temp_file_prefix)s.{region_safe}.tmp.bam | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform=="illumina" and self.caller == "pilon":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "samtools view -T %(ref_file)s -f 0x1 -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && pilon --genome %(ref_file)s --targets {region} --diploid %(calling_params)s --frags %(temp_file_prefix)s.{region_safe}.tmp.bam --variant --output %(temp_file_prefix)s.{region_safe} && bcftools view -e " % vars(self) +  r"""'ALT="."'""" +  " %(temp_file_prefix)s.{region_safe}.vcf | fix_pilon_headers.py --sample %(bam_sample_name)s| add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            elif self.platform=="illumina" and self.caller == "lofreq":
                self.calling_params = calling_params if calling_params else ""
                self.calling_cmd = "samtools view -T %(ref_file)s  -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && lofreq call --call-indels -f %(ref_file)s -r {region} %(calling_params)s  %(temp_file_prefix)s.{region_safe}.tmp.bam  | modify_lofreq_vcf.py --sample %(bam_sample_name)s | add_dummy_AD.py --ref %(ref_file)s --add-dp | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
            else:
                logging.debug("Unknown combination %(platform)s + %(caller)s" % vars(self))
            logging.info("Running variant calling")
            shared_dict['variant_calling'] = self.caller
            run_cmd_parallel_on_genome(self.calling_cmd,ref_file,bed_file = bed_file,threads=threads,desc="Calling variants")
            cmd = "bcftools index  %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self) 
            run_cmd_parallel_on_genome(cmd,ref_file,bed_file = bed_file,threads=threads,desc="Indexing variants")
            temp_vcf_files = ' '.join([f"{tmp}.{r}.vcf.gz" for r in genome_chunks])
            run_cmd(f"bcftools concat -aD {temp_vcf_files} | bcftools view -Oz -o {self.vcf_file}")


        return Vcf(self.vcf_file)
    
    def get_median_depth(
        self,
        ref_file: str,
        software: str = "samtools"
    ):
        logging.info("Calculating median depth using %s" % software)
        shared_dict['depth_calculation'] = software
        if software=="bedtools":
            lines = []
            for l in cmd_out("bedtools genomecov -ibam %s" % (self.bam_file)):
                arr = l.split()
                if arr[0]=="genome":
                    lines.append(arr)
            midpoint =  int(lines[0][3])/2
            x = 0
            for row in lines:
                x = x + int(row[2])
                if x>midpoint:
                    break
            self.median_coverage = int(row[1])
            return int(row[1])
        elif software=="mosdepth":
            self.median_coverage = None
            tmp = str(uuid4())
            run_cmd(f"mosdepth {tmp} {self.bam_file} -f {ref_file}")
            for l in open(f"{tmp}.mosdepth.summary.txt"):
                row = l.strip().split()
                if row[0]=="total":
                    self.median_coverage = int(row[1])
            for f in glob(f"{tmp}*"):
                os.remove(f)
            return int(float(row[3]))
        elif software=="samtools":
            dp = []
            for l in cmd_out(f"samtools depth {self.bam_file}"):
                row = l.strip().split()
                dp.append(int(row[2]))
            ref = FastaFile(ref_file)
            total_len = 0
            for name in ref.references:
                total_len+=ref.get_reference_length(name)
            dp += list([0 for _ in range(total_len - len(dp))])
            return stats.median(dp)

    def calculate_median_coverage(self,ref_file,software="bedtools"):
        logging.info("Calculating median depth")
        if software=="bedtools":
            lines = []
            for l in cmd_out("bedtools genomecov -ibam %s" % (self.bam_file)):
                arr = l.split()

                if arr[0]=="genome":
                    lines.append(arr)
            midpoint =  int(lines[0][3])/2
            x = 0
            for row in lines:
                x = x + int(row[2])
                if x>midpoint:
                    break
            self.median_coverage = int(row[1])
            return int(row[1])
        elif software=="mosdepth":
            self.median_coverage = None
            tmp = str(uuid4())
            run_cmd(f"mosdepth {tmp} {self.bam_file} -f {ref_file}")
            for l in open(f"{tmp}.mosdepth.summary.txt"):
                row = l.strip().split()
                if row[0]=="total":
                    self.median_coverage = int(row[1])
            for f in glob(f"{tmp}*"):
                os.remove(f)
            return int(float(row[3]))


    def get_bed_gt(self,bed_file: str,ref_file: str,caller: str,platform: str):
        logging.info("Getting genotypes for positions in bed file")
        add_arguments_to_self(self, locals())
        results = defaultdict(lambda : defaultdict(dict))
        run_cmd("samtools view -Mb -L %(bed_file)s %(bam_file)s -T %(ref_file)s > %(prefix)s.tmp.bam" % vars(self))
        run_cmd("samtools index %(prefix)s.tmp.bam" % vars(self))
        if platform in ("nanopore","pacbio"):
            caller="bcftools"
        if caller == "gatk":
            cmd = "gatk HaplotypeCaller -I %(prefix)s.tmp.bam -R %(ref_file)s -L %(bed_file)s -OVI false -O /dev/stdout | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        elif caller == "freebayes":
            cmd = "freebayes -f %(ref_file)s -t %(bed_file)s %(prefix)s.tmp.bam --haplotype-length -1 | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        elif caller == "bcftools":
            cmd = "bcftools mpileup -f %(ref_file)s -T %(bed_file)s %(prefix)s.tmp.bam -BI -a AD | bcftools call -mv | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        else:
            cmd = "freebayes -f %(ref_file)s -t %(bed_file)s %(prefix)s.tmp.bam --haplotype-length -1 | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)

        for l in cmd_out(cmd):
            # Chromosome    4348079    0/0    51
            chrom, pos, ref, alt, gt, ad = l.rstrip().split()
            p = GenomePosition(chrom=chrom,pos=int(pos))
            d = {}
            alts = alt.split(",")
            ad = [int(x) for x in ad.split(",")]
            # if gt == "0/0":
            #     d[ref] = ad[0]
            # elif gt == "./.":
            if gt == "./.":
                d[ref] = 0
            else:
                genotypes = list([ref]+alts)
                if platform in ("nanopore","pacbio"):
                    idx = ad.index(max(ad))
                    d[genotypes[idx]] = ad[idx]
                else:
                    for i, a in enumerate(genotypes):
                        d[a] = ad[i]
            results[p] = d

        ref_nt = {}
        for l in cmd_out("bedtools getfasta -fi %s -bed %s" % (ref_file,bed_file)):
            if l[0]==">":
                tmp = l.strip().replace(">","").split(":")
                tmp_chrom = tmp[0]
                tmp_pos = int(tmp[1].split("-")[1])
            else:
                ref_nt[GenomePosition(chrom=tmp_chrom,pos=tmp_pos)] = l.strip().upper()

        for l in cmd_out(f"samtools view -b -L {bed_file} {self.prefix}.tmp.bam | bedtools coverage -a {bed_file} -b - -d -sorted"):
            row = l.strip().split()
            p = GenomePosition(chrom=row[0],pos=int(row[2]))
            cov = int(row[-1])
            if p not in results:
                results[p] = {ref_nt[p]:cov}

        os.remove(f"{self.prefix}.tmp.bam")
        os.remove(f"{self.prefix}.tmp.bam.bai")

        return results

    def calculate_bamstats(self):
        logging.info("Calculating bamstats")
        temp_file = str(uuid4())
        run_cmd(f"samtools stats {self.bam_file} > {temp_file}")
        for l in open(temp_file):
            row = l.strip().split("\t")
            if row[0]!="SN": continue
            if row[1]=='raw total sequences:': self.total_reads = int(row[2])
            if row[1]=='reads mapped:': self.mapped_reads = int(row[2])
        self.pct_reads_mapped = round(self.mapped_reads/self.total_reads*100,2)
        os.remove(temp_file)
    
    def get_missing_genomic_positions(self, bed_file:str, cutoff: int=10) -> List[GenomePositionDepth]:
        """
        Get all positions overlapping bed file that have a depth below a cutoff

        Arguments
        ---------
        bed_file : str
            BED file containing regions to calculate QC metrics for
        cutoff : int
            Depth cutoff to use for calculating QC metrics

        Returns
        -------
        List of GenomePositionDepth objects
        """

        logging.info("Getting missing genomic positions")
        if not hasattr(self,"position_depth"):
            self.calculate_bed_depth(bed_file)
        return [p for p in self.position_depth if p.depth<cutoff]


    def get_region_qc(self,bed_file=None,cutoff=10):
        logging.info("Getting qc metrics for regions")
        if not hasattr(self,"position_depth"):
            self.calculate_bed_depth(bed_file)
        target_qc = []
        for r, data in load_bed(self.bed_file).items():
            target = data[4]
            target_depth = [p for p in self.position_depth if p in r]
            region_len = len(target_depth)
            pos_pass_thresh = len([p for p in target_depth if p.depth>=cutoff])
            target_qc.append(TargetQC(
                target=target, 
                percent_depth_pass=round(pos_pass_thresh/region_len*100,2), 
                median_depth=stats.median([p.depth for p in target_depth]),
            ))
        return target_qc
    
    def get_kmer_counts(self,prefix,klen = 31,threads=1,max_mem=2,counter = "kmc"):
        logging.info("Getting kmer counts")
        if counter=="kmc":
            if threads>32:
                threads = 32
            tmp_prefix = str(uuid4())
            os.mkdir(tmp_prefix)
            bins = "-n128" if platform.system()=="Darwin" else ""
            run_cmd(f"kmc {bins} -fbam -m{max_mem} -t{threads} -k{klen} {self.bam_file} {tmp_prefix} {tmp_prefix}")
            run_cmd(f"kmc_dump {tmp_prefix} {prefix}.kmers.txt")
            run_cmd(f"rm -r {tmp_prefix}*")
            return KmerDump(f"{prefix}.kmers.txt",counter)
        else:
            logging.error("Can't use dsk for bam files, please use kmc instead")

    def get_bam_qc(self, bed_file: str, ref_file:str, depth_cutoff: int, coverage_tool: str = "samtools") -> BamQC:
        """
        Get QC metrics for a bam file
        
        Arguments
        ---------
        bed_file : str
            BED file containing regions to calculate QC metrics for
        ref_file : str
            Reference file to use for calculating median depth
        depth_cutoff : int
            Depth cutoff to use for calculating QC metrics
        coverage_tool : str
            Software to use for calculating median depth. Options are "samtools" or "bedtools"

        Returns
        -------
        BamQC object
        """
        self.calculate_bamstats()
        target_qc = self.get_region_qc(bed_file=bed_file,cutoff=depth_cutoff)
        region_median_depth = stats.median([x.median_depth for x in target_qc])
        genome_median_depth = self.get_median_depth(ref_file=ref_file,software=coverage_tool)
        missing_positions = self.get_missing_genomic_positions(bed_file=bed_file,cutoff=depth_cutoff)
        return BamQC(
            percent_reads_mapped = self.pct_reads_mapped,
            num_reads_mapped = self.mapped_reads,
            target_median_depth = region_median_depth,
            genome_median_depth = genome_median_depth,
            target_qc = target_qc,
            missing_positions = missing_positions,
        )

