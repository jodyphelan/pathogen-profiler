from .utils import iupac
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



def get_barcoding_mutations(mutations: dict, barcode_bed: str) -> List[BarcodePosition]:
    bed = []
    for l in open(barcode_bed):
        if l[0]=="#":
            continue
        row = l.strip().split("\t")
        bed.append(row)
    
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
        snps_report.append(
            BarcodePosition(
                id=marker[3],
                chrom=marker[0],
                pos = marker[2],
                target_allele_count=tmp[1],
                other_allele_count=tmp[0],
                all_allele_count=tmp[0]+tmp[1],
                target_allele_percent=(tmp[1]/sum(tmp))*100
            )
        )
    
        
    return snps_report

def barcode_rows_get_unique_taxa(rows: List[BarcodePosition]) -> List[str]:
    """
    Get unique taxa from a list of BarcodePosition rows.
    
    Arguments
    ---------
    rows: List[BarcodePosition]
        A list of BarcodePosition objects.
    
    Returns
    -------
    List[str]
        A list of unique taxa IDs.
    """
    return sorted(set(row.id for row in rows))

def barcode_rows_calculate_all_allele_count(rows: List[BarcodePosition]) -> None:
    """
    Calculate the total allele count for each BarcodePosition row.
    
    Arguments
    ---------
    rows: List[BarcodePosition]
        A list of BarcodePosition objects.
    
    Returns
    -------
    None
        The function modifies the rows in place.
    """
    for row in rows:
        row.all_allele_count = row.target_allele_count + row.other_allele_count

def barcode_rows_calculate_num_good_sites(rows: List[BarcodePosition],min_percent=2,min_allele_count=5) -> int:
    """
    Calculate the number of good sites based on target allele percent and all allele count.
    A site is considered good if target allele percent >= min_percent and all allele count >= min_allele_count.
    Arguments
    ---------
    rows: List[BarcodePosition]
        A list of BarcodePosition objects.
    Returns
    -------
    int
        The number of good sites.
    """
    num_good_sites = 0
    for row in rows:
        if row.target_allele_percent >= 2 and row.all_allele_count >= 5:
            num_good_sites += 1
    return num_good_sites


def barcode_rows_quantile(rows: List[BarcodePosition], quantile: float) -> float:
    """
    Calculate the quantile of target allele percent for a list of BarcodePosition rows.
    
    Arguments
    ---------
    rows: List[BarcodePosition]
        A list of BarcodePosition objects.
    quantile: float
        The quantile to calculate (between 0 and 1).
    
    Returns
    -------
    float
        The calculated quantile value.
    """
    if not rows:
        return 0.0
    target_allele_percents = [row.target_allele_percent for row in rows]
    return sorted(target_allele_percents)[int(len(target_allele_percents) * quantile)]

def barcode_rows_get_median_frequency(rows: List[BarcodePosition]) -> float:
    """
    Calculate the median frequency of target allele percent for a list of BarcodePosition rows.
    
    Arguments
    ---------
    rows: List[BarcodePosition]
        A list of BarcodePosition objects.
    
    Returns
    -------
    float
        The median frequency of target allele percent.
    """
    if not rows:
        return 0.0
    target_allele_percents = [row.target_allele_percent for row in rows]
    sorted_percents = sorted(target_allele_percents)
    if len(sorted_percents) % 2 == 1:
        return sorted_percents[len(sorted_percents) // 2]
    else:
        mid_index = len(sorted_percents) // 2
        return (sorted_percents[mid_index - 1] + sorted_percents[mid_index]) / 2.0

def barcode(mutations,barcode_bed: str,snps_file=None,iqr_cutoff=15, freq_cutoff=2) -> List[BarcodeResult]:
    bed_num_col = len(open(barcode_bed).readline().rstrip().split("\t"))
    lineage_info = {}
    for l in open(barcode_bed):
        row = l.strip().split("\t")
        lineage_info[row[3]] = row
    
    snps_report = get_barcoding_mutations(mutations,barcode_bed)


    with open(snps_file,"w") if snps_file else open("/dev/null","w") as O:
        for tmp in sorted(snps_report,key=lambda x: x.id):
            O.write("%s\n" % "\t".join([str(x) for x in vars(tmp).values()]))

    barcode_frac = defaultdict(float)
    
    
    barcode_rows_calculate_all_allele_count(snps_report)
    uniq_taxa = barcode_rows_get_unique_taxa(snps_report)
    for taxon in uniq_taxa:
        taxon_rows = [row for row in snps_report if row.id == taxon]
        
        # calculate number of 'good' sites
        num_good_sites = barcode_rows_calculate_num_good_sites(taxon_rows)
        
        if num_good_sites==0:
            logging.debug(f'Skipping {taxon} as no sites pass basic filters')
            continue

        # skip if IQR > cutoff (default=15)
        iqr = barcode_rows_quantile(taxon_rows, 0.75) - barcode_rows_quantile(taxon_rows, 0.25)
        if iqr > iqr_cutoff:
            logging.debug(f'Skipping {taxon} due to high IQR ({iqr})')
            continue

        # skip if median frequency < cutoff (default=2%)
        median_frac = barcode_rows_get_median_frequency(taxon_rows)
        if median_frac < freq_cutoff:
            logging.debug(f'Skipping {taxon} due to low median frequency ({median_frac})')
            continue

        barcode_frac[taxon] = median_frac
        logging.debug(f'Keeping {taxon} with median frequency {median_frac}')

    final_results = []
    for l in barcode_frac:
        tmp = {"id":l,"frequency":barcode_frac[l],"info":[]}
        if bed_num_col>6:
            tmp["info"] = [lineage_info[l][i] for i in range(5,bed_num_col)]
        tmp['support'] = [p for p in snps_report if p.id==tmp['id']]
        final_results.append(BarcodeResult(**tmp))
    return final_results