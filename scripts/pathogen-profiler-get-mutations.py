import pathogenprofiler as pp
import json
import argparse



def main(args):
	if pp.nofolder(args.out_dir):
		pp.run_cmd("mkdir %s" % args.out_dir)
	conf = {
		"ref":args.ref,
		"gff":args.gff,
		"bed":args.bed,
		"ann":args.ann,
		}
	if args.conf:
		conf = json.load(open(args.conf))
	for x in ["ref","gff","bed","ann"]:
		if conf[x]==None:
			pp.log("%s variable is not defined" % x,True)
	bam_obj = pp.bam(
		args.bam,
		args.prefix,
		conf["ref"],
		platform=args.platform
	)
	bcf_obj = bam_obj.call_variants(
		prefix=args.prefix+".targets",
		call_method=args.call_method,
		gff_file=conf["gff"],
		bed_file=conf["bed"],
		mixed_as_missing=False if args.platform == "Illumina" else True,
		threads=args.threads,
		min_dp=args.min_depth,
		af=args.af,
		caller="gatk"
	)
	print("asdhaiduhasuidhui")
	csq = bcf_obj.load_csq(ann_file=conf["ann"])
	variants = []
	chr2gene_pos = {}
	for l in open(conf["ann"]):
		row = l.rstrip().split()
		chr2gene_pos[int(row[1])] = int(row[3])
	for var in list(csq.values())[0]:
		var["_internal_change"] = var["change"]
		var["change"] = pp.reformat_mutations(var["change"],var["type"],var["gene_id"],chr2gene_pos)
		variants.append(var)
	if not args.no_delly:
		delly_bcf = bam_obj.run_delly()
		deletions = delly_bcf.overlap_bed(conf["bed"])
		for deletion in deletions:
			tmp_change = pp.reformat_mutations("%(chr)s_%(start)s_%(end)s" % deletion,var["type"],var["gene_id"],chr2gene_pos)
			tmp = {"genome_pos":deletion["start"],"gene_id":deletion["region"],"chr":deletion["chr"],"freq":1,"type":"large_deletion","change":tmp_change}
			variants.append(tmp)
	json.dump(variants,open("%s/%s.pp-results.json" % (args.out_dir,args.prefix),"w"))
	print("sadjaosidjasiod")
	for x in [".targets.bcf",".targets.csq.bcf",".targets.csq.bcf.csi",".targets.delly.bcf",".targets.delly.bcf.csi",".targets.del_pos.bed",".targets.gvcf.gz",".targets.gvcf.gz.csi",".targets.missing.bcf"]:
		if args.no_delly and "delly" in x: continue
		print(x)
		pp.run_cmd("rm %s%s" % (args.prefix,x))


parser = argparse.ArgumentParser(description='get mutations pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', help='Fasta file')
parser.add_argument('prefix', help='Fasta file')
parser.add_argument('--conf', help='Fasta file')
parser.add_argument('--ref', help='Fasta file')
parser.add_argument('--gff', help='Fasta file')
parser.add_argument('--bed', help='Fasta file')
parser.add_argument('--ann', help='Fasta file')
parser.add_argument('--no-delly',action="store_true", help='Fasta file')
parser.add_argument('--call-method',default="low", help='Fasta file')
parser.add_argument('--platform',default="Illumina", help='Fasta file')
parser.add_argument('--threads',default=1, help='Fasta file')
parser.add_argument('--min-depth',default=10, help='Fasta file')
parser.add_argument('--af',default=0.1, help='Fasta file')
parser.add_argument('--out-dir',default="pp-results", help='Fasta file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
