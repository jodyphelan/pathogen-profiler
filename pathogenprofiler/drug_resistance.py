from collections import defaultdict
from typing import List, Union
from .models import DrVariant, DrGene

def get_lt2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[3]] = row[5].split(",")
    return lt2drugs

def get_gene2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[4]] = row[5].split(",")
    return lt2drugs

def get_drugs2lt(bed_file):
    tmp = get_lt2drugs(bed_file)
    results = defaultdict(list)
    for gene in tmp:
        for drug in tmp[gene]:
            results[drug].append(gene)
    return dict(results)

def get_drugs2gene(bed_file):
    tmp = get_gene2drugs(bed_file)
    results = defaultdict(list)
    for gene in tmp:
        for drug in tmp[gene]:
            results[drug].append(gene)
    return dict(results)

def get_drug_list(bed_file):
    tmp = get_drugs2lt(bed_file)
    return set(tmp.keys())

def get_dr_summary(
    elements: List[Union[DrVariant,DrGene]],
    conf: dict
) -> List[dict]:
    drug_elements = {d:[] for d in conf['drugs']}
    for d in elements:
        for drug in d.get_drugs():
            drug_elements[drug].append(d)

    drug_table = []
    for d in conf['drugs']:
        drug_elements[d] = sorted(drug_elements[d])
        drug_table.append({
            "Drug":d.capitalize(),
            "Genotypic Resistance":"R" if len(drug_elements[d])>0 else "",
            "Mechanisms":", ".join([x.get_str() for x in drug_elements[d]]) if len(drug_elements[d])>0 else ""
        })
    return drug_table
            

def get_summary(json_results,conf,columns = None):
    if not columns:
        columns=[]
    drugs = conf['drugs']
    
    drug_table = []
    results = {}
    annotation = {}

    pool = json_results["dr_variants"].copy()
    if 'resistance_genes' in json_results:
        pool +=  json_results['resistance_genes']
    for x in pool:
        for d in x["drugs"]:
            drug = d["drug"]
            if drug not in results: results[drug] = []
            if "freq" in x:
                results[d["drug"]].append("%s %s (%.2f)" % (x["gene"],x["change"],float(x["freq"])))
            else:
                results[d["drug"]].append("%s (resistance_gene)" % (x["gene"],))
            if drug not in annotation: 
                annotation[drug] = {key:[] for key in columns}
            for key in columns:
                annotation[drug][key].append(d.get(key,""))
    for d in drugs:
        if d in results:
            results[d] = ", ".join(results[d]) if len(results[d])>0 else ""
            r = "R" if len(results[d])>0 else ""
            for key in columns:
                annotation[d][key] = ", ".join(annotation[d][key]) if len(annotation[d][key])>0 else ""
        else:
            results[d] = ""
            r = ""
        dictline = {"Drug":d.capitalize(),"Genotypic Resistance":r,"Mutations":results[d]}
        for key in columns:
            dictline[key] = annotation[d][key] if d in annotation else ""
        drug_table.append(dictline)
    new_json = json_results.copy()
    new_json["drug_table"] = drug_table
    return new_json

def add_drugs_to_variants(variants: List[dict]) -> List[dict]:
    for var in variants:
        if 'annotation' not in var:
            continue
        drug_annotations = [a for a in var['annotation'] if a['type']=='drug_resistance']
        if len(drug_annotations)>0:
            var['drugs'] = drug_annotations
            # for a in var['annotation']:
            #     if a['type']=='drug_resistance':
            #         d = a.copy()
            #         var['drugs'].append(d)
    return variants