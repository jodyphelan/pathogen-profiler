from .utils import iupac
import re
from collections import defaultdict
import logging
from typing import List
from .models import GenomePosition, BarcodeResult, BarcodePosition
import pandas as pd

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


def barcode(mutations,barcode_bed: str,snps_file=None,stdev_cutoff=0.15,iqr=False) -> List[BarcodeResult]:
    if stdev_cutoff is None:
        stdev_cutoff = 0.15
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
    
    rows = []
    for pos in snps_report:
        rows.append(vars(pos))
    df_all = pd.DataFrame(rows)
    df_all.loc[:,'all_allele_count'] = df_all['target_allele_count'] + df_all['other_allele_count']
    
    for taxon in df_all.id.unique():
        df = df_all[df_all.id==taxon].copy()
        pre_filt_num_sites = df.shape[0]
        
        fdf = df.copy() # filtered df

        num_good_sites = df[
            (df['target_allele_percent'] >= 2)
            & (df['all_allele_count'] >= 5)
            
        ].shape[0]
        
        if num_good_sites==0:
            logging.debug(f'Skipping {taxon} as no sites pass basic filters')
            continue

        # skip if number of sites >= 5 and < 25% show alternate
        
        

        # sites_with_alt = fdf[fdf['target_allele_count'] > 0].shape[0]
        # if pre_filt_num_sites>=5 and sites_with_alt/pre_filt_num_sites < 0.25:
        #     logging.debug(f'Skipping {taxon} due to low number of sites ({sites_with_alt}/{pre_filt_num_sites}) with alternate')
        #     continue

        # skip if IQR > 15
        iqr = df['target_allele_percent'].quantile(0.75) - df['target_allele_percent'].quantile(0.25)
        if iqr > 15:
            logging.debug(f'Skipping {taxon} due to high IQR ({iqr})')
            continue

        # skip if median frequency < 2%
        median_frac = df['target_allele_percent'].median()
        if median_frac < 2:
            logging.debug(f'Skipping {taxon} due to low median frequency ({median_frac})')
            continue
        barcode_frac[taxon] = df['target_allele_percent'].median()
        logging.debug(f'Keeping {taxon} with median frequency {median_frac}')

    final_results = []
    for l in barcode_frac:
        tmp = {"id":l,"frequency":barcode_frac[l],"info":[]}
        if bed_num_col>6:
            tmp["info"] = [lineage_info[l][i] for i in range(5,bed_num_col)]
        tmp['support'] = [p for p in snps_report if p.id==tmp['id']]
        final_results.append(BarcodeResult(**tmp))
    return final_results