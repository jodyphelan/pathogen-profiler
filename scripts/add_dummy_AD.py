import sys

for l in sys.stdin:
	if l[0]=="#":
		sys.stdout.write(l)
	else:
		row = l.rstrip().split()
		if "AD" in row[8]:
			sys.stdout.write(l)
		else:
			row[8]+=":AD"
			row[9]+=":0,100"
			sys.stdout.write("%s\n" % "\t".join(row))
