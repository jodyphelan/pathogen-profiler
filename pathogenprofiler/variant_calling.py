import logging
from abc import ABC, abstractmethod
from .utils import cmd_out, run_cmd, run_cmd_parallel_on_genome, load_bed_regions, get_genome_chunks, TempFilePrefix, shared_dict
from uuid import uuid4
from glob import glob
import os
from .vcf import Vcf


class VariantCaller:
    def __init__(
        self, 
        ref_file: str, 
        bam_file: str, 
        prefix: str, 
        bed_file: str = None, 
        threads: int = 1, 
        samclip: bool = False, 
        platform: str = "illumina", 
        calling_params: str = None,
        filters: dict = {}
    ):
        self.temp_file_prefix = str(uuid4())
        self.ref_file = ref_file
        self.bam_file = bam_file
        self.prefix = prefix
        self.bed_file = bed_file
        self.threads = threads
        self.samclip = samclip
        self.platform = platform
        self.calling_params = calling_params if calling_params else ""
        self.af_hard = filters['af_hard'] if 'af_hard' in filters else 0
        if self.platform in ("illumina"):
            self.samclip_cmd = "| samclip --ref %(ref_file)s" % vars(self) if self.samclip else ""
        else:
            self.samclip_cmd = ""
        
        for l in cmd_out("samtools view -H %s" % (bam_file)):
            if l[:3]=="@RG":
                row = l.strip().split("\t")
                for r in row:
                    if r.startswith("SM:"):
                        self.bam_sample_name = r.replace("SM:","")
                        
    @abstractmethod
    def call_variants(self):
        pass

    def run_calling(self, cmd):
        if self.bed_file:
            self.vcf_file = "%s.short_variants.targets.vcf.gz" % (self.prefix) if self.bed_file else "%s.vcf.gz" % (self.prefix)
            genome_chunks = [r.safe for r in load_bed_regions(self.bed_file)]
        else:
            self.vcf_file = "%s.short_variants.vcf.gz" % (self.prefix) if self.bed_file else "%s.vcf.gz" % (self.prefix)
            genome_chunks = [r.safe for r in get_genome_chunks(self.ref_file,self.threads)]

        
        run_cmd_parallel_on_genome(self.calling_cmd,self.ref_file,bed_file = self.bed_file,threads=self.threads,desc="Calling variants")
        cmd = "bcftools index  %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self) 
        run_cmd_parallel_on_genome(cmd,self.ref_file,bed_file = self.bed_file,threads=self.threads,desc="Indexing variants")
        temp_vcf_files = ' '.join([f"{self.temp_file_prefix}.{r}.vcf.gz" for r in genome_chunks])
        run_cmd(f"bcftools concat -aD {temp_vcf_files} | bcftools view -Oz -o {self.vcf_file}")
        for f in glob(self.temp_file_prefix+"*"):
            os.remove(f)

        return Vcf(self.vcf_file)

class BcftoolsCaller(VariantCaller):
    __software__ = "bcftools"
    def call_variants(self) -> Vcf:
        if self.platform=="illumina":
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD,ADF,ADR -r {region} %(temp_file_prefix)s.{region_safe}.tmp.bam | bcftools call -mv | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="nanopore":
            self.calling_cmd = "bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD,ADF,ADR -r {region} %(bam_file)s | bcftools call -mv | annotate_maaf.py | bcftools +fill-tags | bcftools norm -f %(ref_file)s | bcftools filter -e 'IMF < 0.7' -S 0 -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        else:
            raise NotImplementedError("%s not implemented for %s platform" % (self.__software__,self.platform))
        return self.run_calling(self.calling_cmd)

class FreebayesCaller(VariantCaller):
    __software__ = "freebayes"
    def call_variants(self) -> Vcf:
        # Call variants using Freebayes
        if self.platform=="illumina":
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && freebayes -f %(ref_file)s -r {region} --haplotype-length -1 %(calling_params)s %(temp_file_prefix)s.{region_safe}.tmp.bam | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="nanopore":
            self.calling_cmd = "freebayes -f %(ref_file)s -F %(af_hard)s -r {region} --haplotype-length -1 %(calling_params)s %(bam_file)s | annotate_maaf.py  | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        else:
            raise NotImplementedError("%s not implemented for %s platform" % (self.__software__,self.platform))
        return self.run_calling(self.calling_cmd)

class GatkCaller(VariantCaller):
    __software__ = "gatk"
    def call_variants(self) -> Vcf:
        # Call variants using GATK
        if self.platform=="illumina":
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && gatk HaplotypeCaller -R %(ref_file)s -I %(temp_file_prefix)s.{region_safe}.tmp.bam -O /dev/stdout -L {region} %(calling_params)s -OVI false | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        else:
            raise NotImplementedError("%s not implemented for %s platform" % (self.__software__,self.platform))
        return self.run_calling(self.calling_cmd)

class PilonCaller(VariantCaller):
    __software__ = "pilon"
    def call_variants(self) -> Vcf:
        # Call variants using Pilon
        if self.platform=="illumina":
            self.calling_cmd = "samtools view -T %(ref_file)s -f 0x1 -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && pilon --genome %(ref_file)s --targets {region} --diploid %(calling_params)s --frags %(temp_file_prefix)s.{region_safe}.tmp.bam --variant --output %(temp_file_prefix)s.{region_safe} && bcftools view -e " % vars(self) +  r"""'ALT="."'""" +  " %(temp_file_prefix)s.{region_safe}.vcf | fix_pilon_headers.py --sample %(bam_sample_name)s| add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="nanopore":
            self.calling_cmd = """pilon --genome %(ref_file)s --targets {region} %(calling_params)s --nanopore %(bam_file)s --variant --output %(temp_file_prefix)s.{region_safe} && bcftools view -i 'AF>0' %(temp_file_prefix)s.{region_safe}.vcf | fix_pilon_headers.py --sample %(bam_sample_name)s| add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz""" % vars(self)
        else:
            raise NotImplementedError("%s not implemented for %s platform" % (self.__software__,self.platform))
        return self.run_calling(self.calling_cmd)

class LofreqCaller(VariantCaller):
    __software__ = "lofreq"
    def call_variants(self) -> Vcf:
        # Call variants using Lofreq
        if self.platform=="illumina":
            self.calling_cmd = "samtools view -T %(ref_file)s  -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(temp_file_prefix)s.{region_safe}.tmp.bam && samtools index %(temp_file_prefix)s.{region_safe}.tmp.bam && lofreq call --call-indels -f %(ref_file)s -r {region} %(calling_params)s  %(temp_file_prefix)s.{region_safe}.tmp.bam  | modify_lofreq_vcf.py --sample %(bam_sample_name)s | add_dummy_AD.py --ref %(ref_file)s --add-dp | bcftools norm -f %(ref_file)s -Oz -o %(temp_file_prefix)s.{region_safe}.vcf.gz" % vars(self)
        else:
            raise NotImplementedError("%s not implemented for %s platform" % (self.__software__,self.platform))
        return self.run_calling(self.calling_cmd)

