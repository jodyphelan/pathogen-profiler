from __future__ import division
from .bam import *
import os


class fastq:
	"""
	Class for performing mapping to reference.

	Args:
		r1(str): Location of the first read [required]
		r2(str): Location of the second read (if any) use NoneType if there is none.
		ref_file(str): Location of the refrence file [required]
		prefix(str): Prefix for the output files generated (e.g. <prefix>.bam)
		threads(int): Number of threads to use
		platform(str): Platform used in sequencing (currently only 'Illumina' is supported)
		call_method(str): Method used to call variants. Can be ``optimise``, ``low`` or ``high``.
		mapper(str): Mapping tool. Can be ``bwa`` or ``bowtie``

	Returns:
		mapping: A mapping class object
	"""

	def __init__(self,prefix,ref_file,r1,r2=None,threads=4,sample_name=None):
		self.params = {"r1":None,"r2":None}
		self.paired = False
		self.call_method = "high"
		self.mapper = "bwa"
		if r1 and filecheck(r1):
			self.params["r1"] = r1
		else:
			log("Provide at least one fastq file...Exiting\n",True)
		if r2 and filecheck(r2):
			self.params["r2"] = r2
			self.paired = True

		self.params["prefix"] = prefix
		self.params["sample_name"] = sample_name
		self.params["threads"] = threads
		if filecheck(ref_file):
			self.params["ref_file"] = ref_file
		self.params["bam_file"] = "%s.bam" % prefix
		self.params["r1_tp"] = "%s_1_trim_pair.fq" % prefix
		self.params["r1_tu"] = "%s_1_trim_unpair.fq" % prefix
		self.params["r2_tp"] = "%s_2_trim_pair.fq" % prefix
		self.params["r2_tu"] = "%s_2_trim_unpair.fq" % prefix
		self.params["rtu"] = "%s_trim_unpair.fq" % prefix
		self.params["r1t"] = "%s_trim.fq" % prefix

	def trim(self):
		"""Perform trimming"""
		if self.paired:
			cmd = "trimmomatic PE -threads %(threads)s -phred33 %(r1)s %(r2)s %(r1_tp)s %(r1_tu)s %(r2_tp)s %(r2_tu)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % self.params
		else:
			cmd = "trimmomatic SE -threads %(threads)s -phred33 %(r1)s %(r1t)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % self.params
		run_cmd(cmd)
	def illumina(self,mapper="bwa"):
		self.trim()
		if self.paired:
			cmd = "cat %(r1_tu)s %(r2_tu)s > %(rtu)s" % self.params
			run_cmd(cmd)
			psmapper = mapping(self.params["prefix"],self.params["ref_file"],paired1=self.params["r1_tp"],paired2=self.params["r2_tp"],unpaired=self.params["rtu"],mapper=mapper,threads=self.params["threads"],sample_name=self.params["sample_name"])
		else:
			psmapper = mapping(self.params["prefix"],self.params["ref_file"],unpaired=self.params["r1t"],mapper=mapper,threads=self.params["threads"],sample_name=self.params["sample_name"])
		psmapper.map()
		rm_files([self.params["r1t"],self.params["r1_tp"],self.params["r2_tp"],self.params["rtu"],self.params["r1_tu"],self.params["r2_tu"]])
		return psmapper.get_bam(platform="Illumina")
	def minION(self,mapper="minimap2"):
		psmapper = mapping(self.params["prefix"],self.params["ref_file"],unpaired=self.params["r1"],mapper=mapper,threads=self.params["threads"],platform="minION",sample_name=self.params["sample_name"])
		psmapper.map()
		return psmapper.get_bam(platform="minION")
	def get_fastq_qc(self):
		return qc_fastq(self.params["prefix"],self.params["r1"],self.params["r2"])


class mapping:
	"""
	Class for performing mapping to reference.

	Args:
		prefix(str): Prefix for the output files generated (e.g. <prefix>.bam)
		ref_file(str): Location of the refrence file [required]
		paired1(str): First paired end read file
		paired2(str): Seconde paired end read file
		threads(int): Number of threads to use
		platform(str): Platform used in sequencing. Can be ``Illumina`` or ``minION``
		mapper(str): Mapping tool. Can be ``bwa``, ``bowtie`` or ``minimap2``

	Returns:
		mapping: A mapping class object
	"""

	def __init__(self,prefix,ref_file,paired1=None,paired2=None,unpaired=None,threads=4,mapper="bwa",platform="Illumina",sample_name=None):
		if mapper=="bwa": bwa_index(ref_file)
		elif mapper=="bowtie2": bowtie_index(ref_file)
		self.params = {}
		if (paired1 and not paired2) or (paired2 and not paired1):
			log("Please provide two paired end reads...Exiting",True)
		self.paired = True if paired1 and paired2 else False
		if paired1 and filecheck(paired1): self.params["paired1"] = paired1
		if paired2 and filecheck(paired2): self.params["paired2"] = paired2
		if unpaired and filecheck(unpaired): self.params["unpaired"] = unpaired
		self.params["prefix"] = prefix
		self.params["sample_name"] = sample_name if sample_name else prefix
		self.params["threads"] = threads
		self.params["platform"] = platform
		if mapper!="bowtie2" and mapper!="bwa" and mapper!="minimap2":
			log("Please specify correct mapper...Exiting",True)
		self.mapper = mapper
		if filecheck(ref_file): self.params["ref_file"] = ref_file


	def map(self):
		"""Perform mapping"""
		prefix = self.params["prefix"]
		self.params["bwa_prefix"] = "bwa mem -t %(threads)s -c 100 -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -M -T 50" % self.params
		self.params["minimap2_prefix"] = "minimap2 -t %(threads)s -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -a" % self.params
		self.params["bowtie2_prefix"] = "bowtie2 -p %(threads)s --rg-id '%(sample_name)s' --rg 'SM:%(sample_name)s' --rg 'PL:%(platform)s'" % self.params
		self.params["bam_file"] = "%s.bam" % prefix
		if self.params["platform"]=="minION" and self.mapper!="minimap2":
			log("minION data not compatible with %s...Switching to minimap2" % self.mapper)
			self.mapper = "minimap2"
		if self.paired:
			self.params["bam_pair"] = "%s.pair.bam" % prefix
			self.params["bam_single"] = "%s.single.bam" % prefix
			self.params["bam_unsort"] = "%s.unsort.bam" % prefix
			self.params["temp"] = "%s.paired.bam" % self.params["prefix"]
			if self.params["unpaired"]:
				if self.mapper=="bwa":
					cmd = "%(bwa_prefix)s %(ref_file)s %(paired1)s %(paired2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair)s -" % self.params
					run_cmd(cmd)
					cmd = "%(bwa_prefix)s %(ref_file)s %(unpaired)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single)s -" % self.params
					run_cmd(cmd)
				elif self.mapper=="bowtie2":
					cmd = "%(bowtie2_prefix)s -x %(ref_file)s -1 %(paired1)s -2 %(paired2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair)s -" % self.params
					run_cmd(cmd)
					cmd = "%(bowtie2_prefix)s -x %(ref_file)s -U %(unpaired)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single)s" % self.params
					run_cmd(cmd)
				elif self.mapper=="minimap2":
					cmd = "%(minimap2_prefix)s -ax sr %(ref_file)s %(paired1)s %(paired2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair)s -" % self.params
					run_cmd(cmd)
					cmd = "%(minimap2_prefix)s -ax sr %(ref_file)s %(unpaired)s| samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single)s -" % self.params
					run_cmd(cmd)
				cmd = "samtools merge -@ %(threads)s -f %(bam_unsort)s %(bam_pair)s %(bam_single)s" % self.params
				run_cmd(cmd)
				cmd = "samtools sort -@ %(threads)s -o %(bam_file)s %(bam_unsort)s" % self.params
				run_cmd(cmd)
			else:
				if self.mapper=="bwa":
					cmd = "%(bwa_prefix)s %(ref_file)s %(paired1)s %(paired2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
					run_cmd(cmd)
				elif self.mapper=="bowtie2":
					cmd = "%(bowtie2_prefix)s  -x %(ref_file)s -1 %(paired1)s -2 %(paired2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
					run_cmd(cmd)
				elif self.mapper=="minimap2":
					cmd = "%(minimap2_prefix)s -x sr  %(ref_file)s %(paired1)s %(paired2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
					run_cmd(cmd)
			rm_files([self.params["bam_pair"],self.params["bam_single"],self.params["bam_unsort"]])
		else:
			if self.mapper=="bwa":
				cmd = "%(bwa_prefix)s %(ref_file)s %(unpaired)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
			elif self.mapper=="bowtie2":
				cmd = "%(bowtie2_prefix)s  -x %(ref_file)s -U %(unpaired)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s" % self.params
			elif self.mapper=="minimap2":
				if self.params["platform"]=="minION":
					cmd = "%(minimap2_prefix)s -x map-ont  %(ref_file)s %(unpaired)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
				else:
					cmd = "%(minimap2_prefix)s -x sr  %(ref_file)s %(unpaired)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % self.params
			run_cmd(cmd)
	def get_bam(self,platform):
		"""
		Get a bam class object

		Returns:
			bam: A bam class object
		"""
		return bam(self.params["bam_file"],self.params["prefix"],self.params["ref_file"],platform=platform)
