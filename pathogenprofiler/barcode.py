from .utils import infolog, stdev, log, debug, iupac
import json
import re
from collections import defaultdict
from .db import supported_so_terms

def get_missense_codon(x):
    re_obj = re.search("([0-9]+)",x)
    if re_obj:
        return int(re_obj.group(1))
    else:
        log("Error can't find codon number in %s" % x,True)

def get_indel_nucleotide(x):
    re_obj = re.search("([0-9]+)",x)
    if re_obj:
        return int(re_obj.group(1))
    else:
        log("Error can't find nucleotide number in %s" % x,True)



def get_barcoding_mutations(mutations,barcode_bed):
    bed = []
    for l in open(barcode_bed):
        if l[0]=="#":
            continue
        row = l.strip().split("\t")
        bed.append(row)
    
    barcode_support = defaultdict(list)
    snps_report = []
    for marker in bed:
        tmp = [0,0]
        chrom,pos = marker[0],int(marker[2])
        if chrom in mutations and pos in mutations[chrom]:
            for n in iupac(marker[4]):
                if n in mutations[chrom][pos]:
                    tmp[1]+= mutations[chrom][pos][n]
            tmp[0] = sum(list(mutations[chrom][pos].values())) - tmp[1]

        if  tmp==[0,0]: continue
        barcode_support[marker[3]].append(tmp)
        snps_report.append([marker[3],marker[0],marker[2],tmp[1],tmp[0],(tmp[1]/sum(tmp))])
        
    return (barcode_support,snps_report)


def barcode(mutations,barcode_bed,snps_file=None):
    bed_num_col = len(open(barcode_bed).readline().rstrip().split("\t"))
    # bed = []
    lineage_info = {}
    for l in open(barcode_bed):
        row = l.strip().split("\t")
        # bed.append(row)
        lineage_info[row[3]] = row


    barcode_support,snps_report = get_barcoding_mutations(mutations,barcode_bed)
    with open(snps_file,"w") if snps_file else open("/dev/null","w") as O:
        for tmp in sorted(snps_report,key=lambda x: x[0]):
            O.write("%s\n" % "\t".join([str(x) for x in tmp]))

    barcode_frac = defaultdict(float)
    for l in barcode_support:
        # If stdev of fraction across all barcoding positions > 0.15
        # Only look at positions with >5 reads
        tmp_allelic_dp = [x[1]/(x[0]+x[1]) for x in barcode_support[l] if sum(x)>5]
        if len(tmp_allelic_dp)==0: continue
        if stdev(tmp_allelic_dp)>0.15: continue

        # if number of barcoding positions > 5 and only one shows alternate
        if len(barcode_support[l])>5 and len([x for x in barcode_support[l] if (x[1]/(x[0]+x[1]))>0])<2: continue
        barcode_pos_reads = sum([x[1] for x in barcode_support[l]])
        barcode_neg_reads = sum([x[0] for x in barcode_support[l]])
        lf = barcode_pos_reads/(barcode_pos_reads+barcode_neg_reads)
        if lf<0.05:continue
        barcode_frac[l] = lf
    final_results = []

    for l in barcode_frac:
        tmp = {"annotation":l,"freq":barcode_frac[l],"info":[]}
        if bed_num_col>6:
            tmp["info"] = [lineage_info[l][i] for i in range(5,bed_num_col)]
        final_results.append(tmp)
    return final_results

def db_compare(mutations,db):
    annotated_results = mutations
    for i in range(len(mutations["variants"])):
        #var = {'genome_pos': 6140, 'gene_id': 'Rv0005', 'chr': 'Chromosome', 'freq': 0.975609756097561, 'type': 'missense', 'change': '301V>301L'}
        var = mutations["variants"][i]
        for j in range(len(var["consequences"])):
            csq = var["consequences"][j]
            if csq["gene_id"] in db:
                db_var_match = set()
                
                if csq["nucleotide_change"] in db[csq["gene_id"]]:
                    db_var_match.add(json.dumps(db[csq["gene_id"]][csq["nucleotide_change"]]))
                if  csq["protein_change"] in db[csq["gene_id"]]:
                    db_var_match.add(json.dumps(db[csq["gene_id"]][csq["protein_change"]]))
                for t in csq["type"].split("&"):
                    if t in db[csq["gene_id"]]:
                        db_var_match.add(json.dumps(db[csq["gene_id"]][t]))
                if check_for_so_wildcard(csq,db):
                    db_var_match.add(json.dumps(check_for_so_wildcard(csq,db)))
                
                if len(db_var_match)>0:
                    if "annotation" not in annotated_results["variants"][i]["consequences"][j]:
                        annotated_results["variants"][i]["consequences"][j]["annotation"] = []
                    for match in db_var_match:
                        for ann in json.loads(match)["annotations"]:
                            annotated_results["variants"][i]["consequences"][j]["annotation"].append(ann)
    return annotated_results


def extract_affected_positions(change):
    r = re.search("p.[A-Za-z]+(\d+)_[A-Za-z]+(\d+)",change)
    if r:
        return set(range(int(r.group(1)),int(r.group(2))+1))
    r = re.search("p.[A-Za-z]+(\d+)",change)
    if r:
        return set(range(int(r.group(1)),int(r.group(1))+1))
    r = re.search("c.(-*\d+)_(-*\d+)",change)
    if r:
        return set(range(int(r.group(1)),int(r.group(1))+1))
    r = re.search("c.(-*\d+)",change)
    if r:
        return set(range(int(r.group(1)),int(r.group(1))+1))
    r = re.search("n.(-*\d+)_(-*\d+)",change)
    if r:
        return set(range(int(r.group(1)),int(r.group(1))+1))
    r = re.search("n.(-*\d+)",change)
    if r:
        return set(range(int(r.group(1)),int(r.group(1))+1))
    raise ValueError(f"Can't parse {change}")
    
def check_for_so_wildcard(csq,db):
    for var in db[csq['gene_id']]:
        for t in csq['type'].split('&'):
            r = re.search(f"{t}_([pcn])\.(\d+)_(\d+)",var)
            if r: 
                context = r.group(1)
                positions = set(range(int(r.group(2)),int(r.group(3))+1))
                if context=="p":
                    change = csq['protein_change']
                elif context=="c":
                    change = csq['nucleotide_change']
                elif context=="n":
                    change = csq['nucleotide_change']
                affected_positions = extract_affected_positions(change)
                if len(positions.intersection(affected_positions))>0:
                    return db[csq['gene_id']][var]