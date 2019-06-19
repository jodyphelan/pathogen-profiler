import json
from .utils import *
from .bam import *
from .barcode import *
from .fastq import *
from .abi import *
import re
def profiler(conf_file,prefix,r1=None,r2=None,bam_file=None,call_method="low",min_depth=10,platform="Illumina",mapper="bwa",threads=4,run_delly=False,af=0.0,caller="GATK"):
		conf = json.load(open(conf_file))
		for f in conf:
			filecheck(conf[f])
		if not r1 and not r2 and not bam_file:
			log("Please provide at least one fastQ file (-1) or a BAM file (-b)", True)
		elif (r1 or r2) and bam_file:
			log("Please provide fastQ files or a BAM file but not both",True)
		elif not r1 and r2:
			log("Only second fastQ file provided. If profiling a single ended run, just use '-1'",True)
		if r1 or r2:
			if r2:
				fastq_obj = fastq(prefix,conf["ref"],r1,r2,threads=threads,sample_name=prefix.split("/")[-1])
			else:
				fastq_obj = fastq(prefix,conf["ref"],r1,threads=threads,sample_name=prefix.split("/")[-1])
			if platform=="Illumina":
				bam_obj = fastq_obj.illumina(mapper=mapper)
			elif platform=="minION":
				bam_obj = fastq_obj.minION()
		else:
			log("Using %s\n\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)
			bam_obj = bam(bam_file,prefix,conf["ref"],platform=platform)
		bcf_obj = bam_obj.call_variants(prefix=prefix+".targets",call_method=call_method,gff_file=conf["gff"],bed_file=conf["bed"],mixed_as_missing=False if platform == "Illumina" else True,threads=threads,min_dp=min_depth,af=af,caller=caller)
		csq = bcf_obj.load_csq(ann_file=conf["ann"])
		bam_obj.flagstat()


		missing_pos = bcf("%s.targets.missing.bcf" % prefix).get_positions()
		results = {"variants":[],"missing_pos":missing_pos,"qc":{"pct_reads_mapped":bam_obj.pct_reads_mapped,"num_reads_mapped":bam_obj.num_reads_mapped}}
		for sample in csq:
			results["variants"]  = csq[sample]

		mutations = bam_obj.get_bed_gt(conf["barcode"])
		barcode_mutations = barcode(mutations,conf["barcode"])
		results["barcode"] = barcode_mutations
		if run_delly:
			delly_bcf = bam_obj.run_delly()
			deletions = delly_bcf.overlap_bed(conf["bed"])
			for deletion in deletions:
				tmp = {"genome_pos":deletion["start"],"gene_id":deletion["region"],"chr":deletion["chr"],"freq":1,"type":"large_deletion","change":"%(chr)s_%(start)s_%(end)s" % deletion}
				results["variants"].append(tmp)

		results = db_compare(db_file=conf["json_db"],mutations=results)
		delly_pos = []
		for var in results["variants"]:
			if var["type"]=="large_deletion":
				re_obj = re.search("([A-Za-z\.\-\_0-9]+)_([0-9]+)_([0-9]+)",var["change"])
				chrom = re_obj.group(1)
				start = int(re_obj.group(2))
				end = int(re_obj.group(3))
				for i in range(int(start),int(end)+1):
					delly_pos.append((chrom,i))
					if (chrom,i) in results["missing_pos"]:
						results["missing_pos"].remove((chrom,i))
		missing_regions = bam_obj.get_bed_missing(conf["bed"],results["missing_pos"],delly_pos)
		results["missing_regions"] = missing_regions

		return results

def abi_profiler(conf_file,prefix,filenames):
	files = filenames.split(",")
	conf = json.load(open(conf_file))
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
