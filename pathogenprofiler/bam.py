from .utils import add_arguments_to_self, run_cmd, cmd_out, filecheck, index_bam, nofile, rm_files, log, load_bed, median
from .fasta import fasta
from .vcf import bcf, delly_bcf
from tqdm import tqdm
from collections import defaultdict
import sys
class bam:
	"""
	A class to perform operations on BAM files such as SNP calling

	Args:
		bam_file(str): The BAM file [required]
		prefix(str): A prefix for output files [required]
		ref_file(ref_file): A reference (needed by some methods)
		platform(str): Can be either ``Illumina`` or ``minION``
	Returns:
		bam: A bam class object
	"""
	def __init__(self,bam_file,prefix,ref_file,platform="Illumina",threads=4):
		add_arguments_to_self(self,locals())
		index_bam(bam_file,threads=threads)
		filecheck(self.bam_file)
		filecheck(self.ref_file)
		self.ref_fa_obj = fasta(self.ref_file)
		self.ref_fa_dict = self.ref_fa_obj.fa_dict
	def run_delly(self):
		run_cmd("delly call -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self))
		return delly_bcf("%(prefix)s.delly.bcf" % vars(self))
	def call_variants(self,caller="GATK",gff_file=None,bed_file=None,call_method="optimise",min_dp=10,max_dp=None,threads=4,mixed_as_missing=False,low_dp_as_missing=True,af=0.0,whole_genome=False,platform="Illumina",**kwargs):
		add_arguments_to_self(self,locals())
		self.gvcf_file = "%s.gvcf.gz" % (self.prefix if not whole_genome else self.prefix.replace(".targets",""))
		self.missing_vcf_file = "%s.missing.vcf.gz" % self.prefix
		if self.caller=="BCFtools":
			self.bcftools_gbcf(prefix=self.prefix,call_method=call_method,min_dp=min_dp,threads=threads,vtype="both",bed_file=bed_file,low_dp_as_missing=True,platform=platform)
			gbcf_file = self.gvcf_file.replace('.gvcf.gz', '.gbcf')
			self.gbcf_file = gbcf_file
			run_cmd("bcftools view -Oz -l 1 -o %(gvcf_file)s %(gbcf_file)s" % vars(self))
		elif self.caller=="GATK":
			self.gatk_gvcf(bed_file=bed_file if not whole_genome else None, whole_genome=whole_genome)
		else:
			sys.stderr.write("Not sure what to do...")

		self.variant_vcf_file = "%s.vcf.gz" % self.prefix
		self.del_bed = bcf(self.gvcf_file).del_pos2bed(bed_file)
		self.min_dp_cmd = "| bcftools filter -e 'FMT/DP<%(min_dp)s' -Ou -S ." % vars(self) if low_dp_as_missing else ""
		self.max_dp_cmd = "| bcftools filter -e 'FMT/DP>%(max_dp)s' -Ou -S ." % vars(self) if max_dp else ""
		self.mix_cmd = "| bcftools +setGT -- -t q -i 'GT=\"het\" & AD[:1]/(AD[:0]+AD[:1])<0.7' -n . " if mixed_as_missing else ""
		self.af_filter_cmd = "| bcftools view -i 'AF>%s'" % af
		self.bed_option = "-T %s " % bed_file if bed_file else ""


		run_cmd("bcftools view %(bed_option)s %(gvcf_file)s %(min_dp_cmd)s %(max_dp_cmd)s  %(mix_cmd)s %(af_filter_cmd)s | bcftools view -T ^%(del_bed)s -g miss -O z -o %(missing_vcf_file)s" % vars(self))
		run_cmd("bcftools view %(bed_option)s %(gvcf_file)s %(min_dp_cmd)s %(max_dp_cmd)s  %(mix_cmd)s %(af_filter_cmd)s | bcftools view -g ^miss -c 1 -O z -o %(variant_vcf_file)s" % vars(self))

		if gff_file and filecheck(gff_file):
			self.gff_file = gff_file
			self.ann_bcf_file = "%(prefix)s.csq.vcf.gz" % vars(self)
			run_cmd("bcftools csq -p m -f %(ref_file)s -g %(gff_file)s %(variant_vcf_file)s -Oz -o %(ann_bcf_file)s" % vars(self))
			return bcf(self.ann_bcf_file,prefix=self.prefix)
		else:
			return bcf(self.variant_vcf_file,prefix=self.prefix)
	def gatk_gvcf(self,bed_file=None,low_dp_as_missing=False,max_dp=None,min_dp=10,whole_genome=False):
		add_arguments_to_self(self,locals())
		dict_file = self.ref_file.replace(".fasta","").replace(".fa","")+".dict"
		if nofile(dict_file):
			run_cmd("gatk CreateSequenceDictionary -R %(ref_file)s" % vars(self))
		if nofile(self.ref_file+".fai"):
			run_cmd("samtools faidx %(ref_file)s" % vars(self))

		self.gvcf_file = "%s.gvcf.gz" % (self.prefix if not whole_genome else self.prefix.replace(".targets",""))
		self.bed_option = "-L %s " % self.bed_file if self.bed_file else ""
		run_cmd("gatk HaplotypeCaller -R %(ref_file)s -I %(bam_file)s -O %(gvcf_file)s -ERC GVCF %(bed_option)s" % vars(self))
		return bcf(self.gvcf_file)

	def bcftools_gbcf(self,call_method="low",max_dp=None,min_dp=10,threads=4,vtype="snps",bed_file=None,primers=None,overlap_search=True,chunk_size=50000,mpileup_options=None,low_dp_as_missing=False,platform="Illumina",**kwargs):
		"""
		Create a gVCF file (for a description see:https://sites.google.com/site/gvcftools/home/about-gvcf)

		Args:
			ref_file(str): reference file (not required if passed to the bam initiator).
			call_method(str): optimise variant calling based on high or low depth. Options: high|low|optimise
			min_dp(int): Minimum depth required to group site into reference-block
		"""
		add_arguments_to_self(self,locals())
		self.gbcf_file = "%s.gbcf" % self.prefix
		self.cmd_split_chr = "splitchr.py %(ref_file)s %(chunk_size)s --bed %(bed_file)s --reformat" % vars(self) if bed_file else "splitchr.py %(ref_file)s %(chunk_size)s --reformat" % vars(self)
		if primers:
			self.primer_bed_file = "%(prefix)s.primers.bed" % vars(self)
			TMP = open(self.primer_bed_file,"w")
			positions = self.ref_fa.find_primer_positions(primers)
			for x in sorted(positions,key=lambda d:positions[d]["start"]):
				p = positions[x]
				if p["start"] > p["end"]:
					p["start"],p["end"] = p["end"],p["start"]
				TMP.write("%s\t%s\t%s\t%s\n" % (p["chrom"],p["start"],p["end"],x))
			TMP.close()

		if vtype=="snps":
			self.vtype = "-V indels"
		elif vtype=="indels":
			self.vtype = "-V snps"
		elif vtype=="both":
			self.vtype = ""
		else:
			log("Please provide valid vtype: [snps|indels|both]...Exiting!",True)

		self.primer_cmd = " -T ^%(primer_bed_file)s" % vars(self) if primers else ""
		self.extra_cmd = ""
		if call_method=="optimise" and self.platform=="Illumina":
			call_method = self.get_calling_params()
		log("Variant calling optimised for %s" % self.platform)
		self.mpileup_options = ""
		if self.platform=="Illumina" and  call_method=="high":
			self.mpileup_options = "-B -a DP,AD"
		elif self.platform=="Illumina" and call_method=="low":
			self.mpileup_options = "-ABq0 -Q0 -a DP,AD"
		elif self.platform=="minION":
			self.extra_cmd = "|  bcftools filter -e 'IMF < 0.7' -S 0 -Ou"
			if vtype=="snps":
				self.mpileup_options = "-BIq8 -a DP,AD"
			else:
				self.mpileup_options = "-Bq8 -a DP,AD"
		else:
			log("Please choose a valid platform...Exiting!",ext=True)
		if mpileup_options:
			self.mpileup_options = mpileup_options
		self.min_dp_cmd = "| bcftools filter -e 'FMT/DP<%(min_dp)s' -Ou -S ." % vars(self) if low_dp_as_missing else ""
		self.max_dp_cmd = "| bcftools filter -e 'FMT/DP>%(max_dp)s' -Ou -S ." % vars(self) if max_dp else ""
		run_cmd("%(cmd_split_chr)s | parallel --col-sep '\\t' -j %(threads)s \"bcftools mpileup  -f %(ref_file)s %(bam_file)s %(mpileup_options)s -r {1} | bcftools call %(primer_cmd)s %(vtype)s -mg %(min_dp)s | bcftools norm -f %(ref_file)s %(min_dp_cmd)s %(max_dp_cmd)s %(extra_cmd)s | bcftools view -Ob -o %(prefix)s_{2}.bcf \"" % vars(self))
		run_cmd("%(cmd_split_chr)s | awk '{print \"%(prefix)s_\"$2\".bcf\"}' | parallel -j  %(threads)s \"bcftools index {}\"" % vars(self))

		if primers:
			self.non_primer_bcf = "%(prefix)s.non_primer.bcf" % vars(self)
			self.primer_bcf = "%(prefix)s.primer.bcf" % vars(self)
			if overlap_search:
				self.generate_primer_bcf()
			else:
				run_cmd("bcftools mpileup  -f %(ref_file)s %(bam_file)s %(mpileup_options)s -B -R %(primer_bed_file)s | bcftools call %(vtype)s -m | bcftools +setGT -Ob -o %(primer_bcf)s -- -t a -n ." % vars(self))
			run_cmd("bcftools concat -aD -Ob -o %(non_primer_bcf)s `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$2\".bcf\"}'`" % vars(self))
			run_cmd("bcftools concat %(primer_bcf)s %(non_primer_bcf)s | bcftools sort -Ob -o %(gbcf_file)s " % vars(self))
		else:
			run_cmd("bcftools concat -aD -Ob -o %(gbcf_file)s `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$2\".bcf\"}'`" % vars(self))

		run_cmd("rm `%(cmd_split_chr)s  | awk '{print \"%(prefix)s_\"$2\".bcf*\"}'`" % vars(self))
		if primers:
			rm_files([self.non_primer_bcf,self.primer_bcf])
		return bcf(self.gbcf_file,prefix=self.prefix)
	def get_calling_params(self):
		dp = []
		log("Optimising call method")
		for l in cmd_out("samtools depth %(bam_file)s" % vars(self)):
			arr = l.rstrip().split()
			dp.append(int(arr[2]))
		med_dp = median(dp)
		log("Median depth: %s" % med_dp)
		if med_dp<30:
			log("Using low depth approach")
			return "low"
		else:
			log("Using high depth approach")
			return "high"
	def flagstat(self):
		lines = []
		for l in cmd_out("samtools flagstat %s" % (self.bam_file)):
			arr = l.split()
			lines.append(arr)
		self.num_reads_mapped = int(lines[4][0])
		self.pct_reads_mapped = 0.0 if self.num_reads_mapped==0 else float(lines[4][4][1:-1])
		return self.num_reads_mapped,self.pct_reads_mapped
	def load_genome_cov(self,bed_file=None):
		add_arguments_to_self(self,locals())
		self.ref_cov = {}
		for s in self.ref_fa_dict:
			self.ref_cov[s] = [0 for x in range(len(self.ref_fa_dict[s]))]
		if bed_file:
			samtools_cmd = "samtools view -bL %s %s | samtools depth -aa -b %s --reference %s -" % (bed_file,self.bam_file,bed_file,self.ref_file)
		else:
			samtools_cmd = "samtools depth -aa --reference %s %s" % (self.ref_file,self.bam_file)
		for line in cmd_out(samtools_cmd):
			arr = line.split()
			if arr[0] not in self.ref_cov: log("Can't find %s in FASTA...Have you used the correct reference sequence?" % arr[0]);quit()
			self.ref_cov[arr[0]][int(arr[1])-1] = int(arr[2])
	def get_bed_missing(self,bed_file,missing_pos=None,deletions_pos=None):
		# BED lines like this: Chromosome	10	20	region1
		bed_pos = load_bed(bed_file,[1,2,3,4],4,intasint=True)
		self.load_genome_cov(bed_file)
		miss_region = {}
		for region in bed_pos:
			miss_pos = 0
			chrom,start,end = bed_pos[region][0:3]
			start = int(start)
			end = int(end)
			for i in range(start, end):
				if ((chrom,i) in missing_pos or self.ref_cov[chrom][i]==0) and (chrom,i) not in deletions_pos:
					miss_pos+=1
			miss_region[region] = miss_pos/(end-start)
		return miss_region
	def get_bed_gt(self,bed_file,caller="GATK",platform="Illumina"):
		add_arguments_to_self(self, locals())
		results = defaultdict(lambda : defaultdict(dict))
		if caller=="GATK":
			cmd = "gatk HaplotypeCaller -I %(bam_file)s -R %(ref_file)s -L %(bed_file)s -ERC BP_RESOLUTION -OVI false -O /dev/stdout | bcftools view -a | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
		else:
			cmd = "bcftools mpileup -f %(ref_file)s -R %(bed_file)s %(bam_file)s -BI -a AD | bcftools call -m | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
		for l in cmd_out(cmd):
			#Chromosome	4348079	0/0	51
			chrom,pos,ref,alt,gt,ad = l.rstrip().split()
			pos =int(pos)
			d = {}
			alts = alt.split(",")
			ad = [int(x) for x in ad.split(",")]
			if gt=="0/0":
				d[ref] = ad[0]
			elif gt=="./.":
				d[ref] = 0
			else:
				genotypes = list([ref]+alts)
				if platform=="Illumina":
					idx = genotypes.index(max(genotypes))
					d[genotypes[idx]] = ad[idx]
				else:
					for i,a in enumerate(genotypes):
						d[a] = ad[i]
			results[chrom][pos] = d
		return results
