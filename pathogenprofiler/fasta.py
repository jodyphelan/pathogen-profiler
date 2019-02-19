from collections import OrderedDict
class fasta:
	"""
	Class to represent fasta seuqnces in a python dict.

	Args:
		filename(str): Location of the fasta file

	Returns:
		fasta: A fasta class object
	"""
	def __init__(self,filename):
		fa_dict = OrderedDict()
		seq_name = ""
		self.fa_file = filename
		for l in open(filename):
			line = l.rstrip()
			if line.startswith(">"):
				seq_name = line[1:].split()[0]
				fa_dict[seq_name] = []
			else:
				fa_dict[seq_name].append(line)
		result = {}
		counter = 0
		sum_length = {}
		for seq in fa_dict:
			result[seq] = "".join(fa_dict[seq])
			result[seq] = result[seq].upper()
			sum_length[(counter+1,counter+len(result[seq]))] = seq
			counter = counter+len(result[seq])
		self.sum_length = sum_length
		self.fa_dict = result
