from .utils import filecheck, log, run_cmd
from .bam import bam
from .barcode import barcode, db_compare
from .fastq import fastq
from .abi import abi
from .vcf import bcf
from .fasta import fasta
import re


def profiler(conf, prefix, r1=None, r2=None, bam_file=None, call_method="low", min_depth=10, platform="Illumina", mapper="bwa", threads=4, run_delly=False, af=0.0, caller="GATK", whole_genome=False,notrim=False,noflagstat=False):
		for f in conf:
			filecheck(conf[f])
		if not r1 and not r2 and not bam_file:
			log("Please provide at least one fastQ file (-1) or a BAM file (-b)", True)
		elif (r1 or r2) and bam_file:
			log("Please provide fastQ files or a BAM file but not both", True)
		elif not r1 and r2:
			log("Only second fastQ file provided. If profiling a single ended run, just use '-1'",True)
		if r1 or r2:
			if r2:
				fastq_obj = fastq(prefix, conf["ref"], r1, r2, threads=threads, sample_name=prefix.split("/")[-1])
			else:
				fastq_obj = fastq(prefix, conf["ref"], r1, threads=threads, sample_name=prefix.split("/")[-1])
			if platform == "Illumina":
				bam_obj = fastq_obj.illumina(mapper=mapper,notrim=notrim)
			elif platform == "minION":
				bam_obj = fastq_obj.minION()
		else:
			log("Using %s\n\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)
			bam_obj = bam(bam_file, prefix, conf["ref"], platform=platform)
		bcf_obj = bam_obj.call_variants(prefix=prefix+".targets", whole_genome=whole_genome,call_method=call_method, gff_file=conf["gff"], bed_file=conf["bed"], mixed_as_missing=False if platform == "Illumina" else True, threads=threads, min_dp=min_depth, af=af, caller=caller if platform =="Illumina" else "BCFtools", platform=platform)
		csq = bcf_obj.load_csq(ann_file=conf["ann"])
		if noflagstat:
			bam_obj.pct_reads_mapped = "NA"
			bam_obj.num_reads_mapped = "NA"
		else:
			bam_obj.flagstat()


		missing_pos = bcf("%s.targets.missing.vcf.gz" % prefix).get_positions()
		results = {"variants":[],"missing_pos":missing_pos,"qc":{"pct_reads_mapped":bam_obj.pct_reads_mapped,"num_reads_mapped":bam_obj.num_reads_mapped}}
		for sample in csq:
			results["variants"]  = csq[sample]
		mutations = bam_obj.get_bed_gt(conf["barcode"],caller="BCFtools" if (caller=="BCFtools" or platform=="minION") else "BCFtools")
		barcode_mutations = barcode(mutations,conf["barcode"])
		results["barcode"] = barcode_mutations
		if platform=="minION":
			run_delly = False
		if run_delly:
			delly_bcf = bam_obj.run_delly()
			deletions = delly_bcf.overlap_bed(conf["bed"])
			for deletion in deletions:
				tmp = {
					"genome_pos": deletion["start"], "gene_id": deletion["region"],
					"chr": deletion["chr"], "freq": 1, "type": "large_deletion",
					"change": "%(chr)s_%(start)s_%(end)s" % deletion
					}
				results["variants"].append(tmp)

		results = db_compare(db_file=conf["json_db"], mutations=results)
		delly_pos = []
		for var in results["variants"]:
			if var["type"] == "large_deletion":
				re_obj = re.search("([A-Za-z\.\-\_0-9]+)_([0-9]+)_([0-9]+)", var["change"])
				chrom = re_obj.group(1)
				start = int(re_obj.group(2))
				end = int(re_obj.group(3))
				for i in range(int(start),int(end)+1):
					delly_pos.append((chrom,i))
					if (chrom, i) in results["missing_pos"]:
						results["missing_pos"].remove((chrom,i))
		missing_regions = bam_obj.get_bed_missing(conf["bed"], results["missing_pos"], delly_pos)
		results["missing_regions"] = missing_regions

		return results


def fasta_profiler(conf, prefix, filename):
	fasta_obj = fasta(filename)
	wg_vcf_file = fasta_obj.get_ref_variants(conf["ref"], prefix, gff=conf["gff"])
	wg_vcf_obj = bcf(wg_vcf_file)
	targets_vcf_file = prefix+".targets.vcf.gz"
	run_cmd("bcftools view -c 1 %s -Oz -o %s -T %s" % (wg_vcf_file,targets_vcf_file,conf['bed']))
	targets_vcf_obj = bcf(targets_vcf_file)
	csq = targets_vcf_obj.load_csq(ann_file=conf["ann"])
	results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
	for sample in csq:
		results["variants"]  = csq[sample]
	mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
	barcode_mutations = barcode(mutations,conf["barcode"])
	results["barcode"] = barcode_mutations
	results = db_compare(db_file=conf["json_db"], mutations=results)
	return results


def abi_profiler(conf,prefix,filenames):
	files = filenames.split(",")
	for f in conf:
		filecheck(conf[f])
	for f in files:
		filecheck(f)
	abi_obj = abi(files,prefix)
	vcf_obj = abi_obj.get_variants_vcf(conf["ref"],conf["gff"])
	csq = vcf_obj.load_csq(ann_file=conf["ann"])
	results = {"variants":[]}
	for sample in csq:
		results["variants"]  = csq[sample]
	return results
