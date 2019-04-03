import sys
import gzip
import os.path
import subprocess
import csv
from collections import defaultdict
import json
import random
import math
rand_generator = random.SystemRandom()

def filetype(x):
	for l in cmd_out("file %s" % x):
		pass
	row = l.rstrip().split()
	return row[1]

def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])

def stdev(arr):
	mean = sum(arr)/len(arr)
	return math.sqrt(sum([(x-mean)**2 for x in arr])/len(arr))

def add_arguments_to_self(self,args):
	for x in args:
		if x=="self": continue
		vars(self)[x] = args[x]
	if "kwargs" in args:
		for x in args["kwargs"]:
			vars(self)[x] = args["kwargs"][x]

def cmd_out(cmd,verbose=1):
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stderr = open("/dev/null","w")
	else:
		stderr = open("/dev/null","w")
	try:
		res = subprocess.Popen(cmd,shell=True,stderr = stderr,stdout=subprocess.PIPE)
		for l in res.stdout:
			yield l.decode().rstrip()
	except:
		print("Command Failed! Please Check!")
		exit(1)
	stderr.close()

def get_random_file(prefix = None,extension=None):
	randint = rand_generator.randint(1,999999)
	if prefix:
		if extension:
			return "%s.%s.%s" % (prefix,randint,extension)
		else:
			return "%s.%s.txt" % (prefix,randint)
	else:
		if extension:
			return "%s.tmp.%s" % (randint,extension)
		else:
			return "%s.tmp.txt" % (randint)

def log(msg,ext=False):
	sys.stderr.write("\n"+str(msg)+"\n")
	if ext:
		exit(1)

def init_params():
	conf = json.load(open("%s/%s" % (sys.prefix,"pathogenseq.conf")))
	return conf

def load_tsv(filename):
	meta = {}
	for row in csv.DictReader(open(filename),delimiter="\t"):
		if "sample" not in row:
			print("No sample column...Exiting")
			quit(1)
		meta[row["sample"]] = {}
		columns = set(row)-set(["sample"])
		for c in columns:
			meta[row["sample"]][c.upper()] = row[c]
	columns = [c.upper() for c in set(row)-set(["sample"])]
	return columns,meta

def load_bed(filename,columns,key1,key2=None,intasint=False):
	results = defaultdict(lambda: defaultdict(tuple))
	for l in open(filename):
		row = l.rstrip().split("\t")
		if key2:
			if max(columns+[key1,key2])>len(row):
				log("Can't access a column in BED file. The largest column specified is too big",True)
			if key2==2 or key2==3:
				results[row[key1-1]][int(row[key2-1])] = tuple([row[int(x)-1] for x in columns])
			else:
				results[row[key1-1]][row[key2-1]] = tuple([row[int(x)-1] for x in columns])
		else:
			if max(columns+[key1])>len(row):
				log("Can't access a column in BED file. The largest column specified is too big",True)
			results[row[key1-1]]= tuple([row[int(x)-1] for x in columns])
	return results

def split_bed(bed_file,size,reformat=False):
	bed_regions = {}
	for l in open(filecheck(bed_file)):
		row = l.rstrip().split()
		chrom,start,end = row[0],int(row[1]),int(row[2])
		if end-start>size:
			tmps = start
			tmpe = start
			while tmpe<end:
				tmpe+=size
				if tmpe>end:
					tmpe=end
				loc = "%s:%s-%s" % (chrom,tmps,tmpe)
				loc_str = "%s_%s_%s" % (chrom,tmps,tmpe)
				if reformat:
					sys.stdout.write("%s\t%s\n" % (loc,loc_str))
				else:
					sys.stdout.write("%s\n"%loc)
				tmps=tmpe+1
		else:
			loc = "%s:%s-%s" % (chrom,start,end)
			loc_str = "%s_%s_%s" % (chrom,start,end)
			if reformat:
				sys.stdout.write("%s\t%s\n" % (loc,loc_str))
			else:
				sys.stdout.write("%s\n"%loc)
def filecheck(filename):
	"""
	Check if file is there and quit if it isn't
	"""
	if not os.path.isfile(filename):
		print("Can't find %s" % filename)
		exit(1)
	else:
		return filename

def foldercheck(filename):
	"""
	Check if file is there and quit if it isn't
	"""
	if not os.path.isdir(filename):
		print("Can't find %s" % filename)
		exit(1)
	else:
		return filename

def debug(s):
	sys.stderr.write("#"*40+"\n")
	sys.stderr.write("%s\n" % s)
	sys.stderr.write("#"*40+"\n")

def nofile(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isfile(filename):
		return True
	else:
		return False

def nofolder(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isdir(filename):
		return True
	else:
		return False

def bowtie_index(ref):
	if nofile("%s.1.bt2"%ref):
		cmd = "bowtie2-build %s %s" % (ref,ref)
		run_cmd(cmd)

def bwa_index(ref):
	"""
	Create BWA index for a reference
	"""
	if nofile("%s.bwt"%ref):
		cmd = "bwa index %s" % ref
		run_cmd(cmd)

def run_cmd(cmd,verbose=1,target=None):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	if target and filecheck(target): return True
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/stdout","w")
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")
	else:
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		print("Command Failed! Please Check!")
		exit(1)

def index_bam(bamfile,threads=4,overwrite=False):
	"""
	Indexing a bam file
	"""
	cmd = "samtools index -@ %s %s" % (threads,bamfile)
	if filecheck(bamfile):
		if nofile(bamfile+".bai"):
			run_cmd(cmd)
		elif os.path.getmtime(bamfile+".bai")<os.path.getmtime(bamfile) or overwrite:
			run_cmd(cmd)

def index_bcf(bcffile,threads=4,overwrite=False):
	"""
	Indexing a bam file
	"""
	cmd = "bcftools index --threads %s -f %s" % (threads,bcffile)
	if filecheck(bcffile):
		if nofile(bcffile+".csi"):
			run_cmd(cmd)
		elif os.path.getmtime(bcffile+".csi")<os.path.getmtime(bcffile) or overwrite:
			run_cmd(cmd)

def verify_fq(filename):
	"""
	Return True if input is a valid fastQ file
	"""
	FQ = open(filename) if filename[-3:]!=".gz" else gzip.open(filename)
	l1 = FQ.readline()
	if l1[0]!="@":
		print("First character is not \"@\"\nPlease make sure this is fastq format\nExiting...")
		exit(1)
	else:
		return True

def rm_files(x,verbose=True):
	"""
	Remove a files in a list format
	"""
	for f in x:
		if os.path.isfile(f):
			if verbose: print("Removing %s" % f)
			os.remove(f)

def file_len(filename):
	"""
	Return length of a file
	"""
	filecheck(filename)
	for l in subprocess.Popen("wc -l %s" % filename,shell=True,stdout=subprocess.PIPE).stdout:
		res = l.rstrip().split()[0]
	return int(res)

def gz_file_len(filename):
	"""
	Return lengths of a gzipped file
	"""
	filecheck(filename)
	for l in subprocess.Popen("gunzip -c %s |wc -l" % filename,shell=True,stdout=subprocess.PIPE).stdout:
		res = l.rstrip().split()[0]
	return int(res)

def download_from_ena(acc):
	if len(acc)==9:
		dir1 = acc[:6]
		cmd = "wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s*" % (dir1,acc,acc)
	elif len(acc)==10:
		dir1 = acc[:6]
		dir2 = "00"+acc[-1]
		cmd = "wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s*" % (dir1,dir2,acc,acc)
	else:
		print("Check Accession: %s" % acc)
		exit(1)
	run_cmd(cmd)

def which(program):
	import os
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

def programs_check(programs):
	for p in programs:
		if which(p)==None:
			log("Can't find %s in path... Exiting." % p)
			quit(1)
