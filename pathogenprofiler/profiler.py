import json
from .utils import *
from .bam import *
from .barcode import *
from .fastq import *

def profiler(conf_file,prefix,r1=None,r2=None,bam_file=None,call_method="low",min_depth=10,platform="Illumina",mapper="bwa",threads=4,run_delly=False):
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
		bcf_obj = bam_obj.call_variants(prefix=prefix+".targets",call_method=call_method,gff_file=conf["gff"],bed_file=conf["bed"],mixed_as_missing=False if platform == "Illumina" else True,threads=threads,min_dp=min_depth)
		csq = bcf_obj.load_csq(ann_file=conf["ann"])
		bam_obj.flagstat()


		missing_pos = bcf("%s.targets.missing.bcf" % prefix).get_positions()
		missing_regions = bam_obj.get_bed_missing(conf["bed"],missing_pos)
		results = {"variants":[],"missing_regions":missing_regions,"missing_pos":missing_pos,"qc":{"pct_reads_mapped":bam_obj.pct_reads_mapped,"num_reads_mapped":bam_obj.num_reads_mapped}}
		for sample in csq:
			results["variants"]  = csq[sample]

		mutations = bam_obj.get_bed_gt(conf["barcode"])
		barcode_mutations = barcode(mutations,conf["barcode"])
		results["barcode"] = barcode_mutations
		results = db_compare(db_file=conf["json_db"],mutations=results)
		if run_delly:
			delly_bcf = bam_obj.run_delly()
			deletions = delly_bcf.overlap_bed(conf["bed"])
			results = annotate_deletions(deletions = deletions, mutations=results,bed_file=conf["bed"])
		return results
