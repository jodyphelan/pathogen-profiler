from .utils import stdev, iupac
import re
from collections import defaultdict
import logging
from typing import List
from .models import GenomePosition, BarcodeResult, BarcodePosition

def get_missense_codon(x):
    re_obj = re.search("([0-9]+)",x)
    if re_obj:
        return int(re_obj.group(1))
    else:
        logging.error("Error can't find codon number in %s" % x,True)

def get_indel_nucleotide(x):
    re_obj = re.search("([0-9]+)",x)
    if re_obj:
        return int(re_obj.group(1))
    else:
        logging.error("Error can't find nucleotide number in %s" % x,True)



def get_barcoding_mutations(mutations: dict, barcode_bed: str) -> tuple[dict, List[list]]:
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
        p = GenomePosition(chrom=marker[0],pos=int(marker[2]))
        if p in mutations:
            for n in iupac(marker[4]):
                if n in mutations[p]:
                    tmp[1]+= mutations[p][n]
            tmp[0] = sum(list(mutations[p].values())) - tmp[1]

        if  tmp==[0,0]: continue
        barcode_support[marker[3]].append(tmp)
        snps_report.append(
            BarcodePosition(
                id=marker[3],
                chrom=marker[0],
                pos = marker[2],
                target_allele_count=tmp[1],
                other_allele_count=tmp[0],
                target_allele_percent=(tmp[1]/sum(tmp))*100
            )
        )
    
        
    return (barcode_support,snps_report)


def barcode(mutations,barcode_bed: str,snps_file=None,stdev_cutoff=0.15) -> List[BarcodeResult]:
    bed_num_col = len(open(barcode_bed).readline().rstrip().split("\t"))
    # bed = []
    lineage_info = {}
    for l in open(barcode_bed):
        row = l.strip().split("\t")
        # bed.append(row)
        lineage_info[row[3]] = row


    barcode_support,snps_report = get_barcoding_mutations(mutations,barcode_bed)

    with open(snps_file,"w") if snps_file else open("/dev/null","w") as O:
        for tmp in sorted(snps_report,key=lambda x: x.id):
            O.write("%s\n" % "\t".join([str(x) for x in vars(tmp).values()]))

    barcode_frac = defaultdict(float)
    for l in barcode_support:
        logging.debug("Processing %s" % l)
        logging.debug(barcode_support[l])
        # If stdev of fraction across all barcoding positions > stdev_cutoff
        # Only look at positions with >5 reads
        tmp_allelic_dp = [x[1]/(x[0]+x[1]) for x in barcode_support[l] if sum(x)>5]
        logging.debug(tmp_allelic_dp)

        # remove positions with no SNP (fraction=0)
        tmp_allelic_dp = [x for x in tmp_allelic_dp if x>0]
        logging.debug(tmp_allelic_dp)

        if len(tmp_allelic_dp)==0: 
            logging.debug("No SNPs found for %s" % l)
            continue
        if stdev(tmp_allelic_dp)>stdev_cutoff:
            logging.debug(f"Stdev {stdev(tmp_allelic_dp)} > {stdev_cutoff} for {l}")
            continue

        # if number of barcoding positions > 5 and only one shows alternate
        num_positions_with_alt = len([x for x in barcode_support[l] if (x[1]/(x[0]+x[1]))>0])
        logging.debug("Number of positions with alternate %s" % num_positions_with_alt)
        if len(barcode_support[l])>5 and num_positions_with_alt<2:
            logging.debug("Number of positions with alternate <2 for %s" % l)
            continue
        # if number of barcoding positions > 5 and there are less than 5% of possible positions with alternate
        if len(barcode_support[l])>5 and num_positions_with_alt<=0.10*len(barcode_support[l]):
            logging.debug("Number of positions with alternate < 10 percent for %s" % l)
            continue
        
        barcode_pos_reads = sum([x[1] for x in barcode_support[l]])
        barcode_neg_reads = sum([x[0] for x in barcode_support[l]])
        lf = barcode_pos_reads/(barcode_pos_reads+barcode_neg_reads)
        logging.debug("Barcode %s has frequency %s" % (l,lf))
        if lf<0.05:
            logging.debug("Frequency < 0.05 for %s" % l)
            continue
        barcode_frac[l] = lf
    final_results = []

    for l in barcode_frac:
        tmp = {"id":l,"frequency":barcode_frac[l],"info":[]}
        if bed_num_col>6:
            tmp["info"] = [lineage_info[l][i] for i in range(5,bed_num_col)]
        tmp['support'] = [p for p in snps_report if p.id==tmp['id']]
        final_results.append(BarcodeResult(**tmp))
    return final_results