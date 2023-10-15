from collections import defaultdict

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