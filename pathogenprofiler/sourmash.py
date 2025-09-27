
import logging
from typing import List
from .utils import run_cmd
import csv
from uuid import uuid4
import os
from .models import TaxonomicHit, Species


def combine_species_abundance(matches:List[dict]) -> List[dict]:
    species_detected = set(t.species for t in matches)
    species_objects = []
    if len(matches) == 0:
        return []

    for species in species_detected:
        hits = [t for t in matches if t.species == species]
        hits = sorted(hits,key=lambda x: x.abundance,reverse=True)
        species_objects.append(hits[0])

    total_abundance = sum([s.abundance for s in species_objects])
    for s in species_objects:
        s.relative_abundance = s.abundance/total_abundance*100

    species_objects = sorted(species_objects,key=lambda x: x.relative_abundance,reverse=True)
    return species_objects

class Sketch:

    def __init__(self,filename, tmp_prefix=None):
        self.filename = filename
        self.tmp_prefix = tmp_prefix if tmp_prefix else str(uuid4())

    def get_species_hits(self, ref_db, db_annotation, *args, **kwargs):
        sequence_hits = self.classify(ref_db, db_annotation, *args, **kwargs)
        print("adjoaidjadoiajd")
        accession_data = {}
        for row in csv.DictReader(open(db_annotation)):
            accession_data[row["accession"]] = row

        species_hits = []
        for hit in sequence_hits:
            data = hit.model_dump()
            data.update(accession_data.get(hit.accession,{}))
            species_hits.append(Species(**data))
            # if hit.accession in accession_data:
            #     hit.species = accession_data[hit.accession]["species"]
            #     hit.ncbi_organism_name = accession_data[hit.accession]["ncbi_organism_name"]

        combined_species_hits = combine_species_abundance(species_hits)
        print(combined_species_hits)
        quit()
    


class SourmashSig(Sketch):

    def __init__(self,filename,tmp_prefix=None):
        super().__init__(filename, tmp_prefix=tmp_prefix)

    # def filter(self,min_abundance=5,outfile=None):
    #     logging.info("Filtering sourmash sig")
        
    #     if outfile is None:
    #         outfile = self.filename.replace(".sig",".filtered.sig")
        
    #     run_cmd(f"sourmash sig filter -o {outfile} {self.filename} -m {min_abundance}")
        
    #     return SourmashSig(outfile,tmp_prefix=self.tmp_prefix)
    
    # def search(self, ref_db, db_annotation, ani_threshold=95):
    #     logging.info("Searching sourmash sig")
        
    #     outfile = "%s" % self.tmp_prefix+".sourmash.csv"
    #     run_cmd(f"sourmash search -n 10 {self.filename} {ref_db} --containment --estimate-ani --ignore-abundance -o {outfile}")

    #     accession_data = {}
    #     for row in csv.DictReader(open(db_annotation)):
    #         accession_data[row["accession"]] = row

    #     results = []
    #     for row in csv.DictReader(open(outfile)):
    #         d = {
    #             'accession': row['name'],
    #             'ani': float(row['ani']),
    #         }
    #         d.update(accession_data[row['name']])
    #         results.append(d)
    #     results = [x for x in results if x["ani"]>=ani_threshold]
    #     return results[:10]

    def classify(self, ref_db, db_annotation=None, intersect_bp=500000,f_match_threshold=0.1):
        logging.info("Classifying sourmash sig")

        outfile = "%s" % self.tmp_prefix+".sourmash.csv"
        run_cmd(f"sourmash gather {self.filename} {ref_db} -o {outfile}")

       

        results = []
        filtered_rows = []
        taxonomic_hits = []
        if not os.path.exists(outfile):
            return []
        for row in csv.DictReader(open(outfile)):
            logging.debug(row)
            if intersect_bp and float(row['intersect_bp'])<intersect_bp:
                logging.debug("skipping because of intersect_bp")
                continue
            if f_match_threshold and float(row['f_match'])<f_match_threshold:
                logging.debug("skipping because of f_match")
                continue

            filtered_rows.append(row)
            hit = TaxonomicHit(
                prediction_method="sourmash",
                accession=row["name"],
                ani=round(float(row["match_containment_ani"])*100,2),
                abundance=float(row["average_abund"])
            )
            taxonomic_hits.append(hit)

        return taxonomic_hits
    


class SylphSketch:
    def classify(self, ref_db, *args, **kwargs):
        outfile = "%s" % self.tmp_prefix+".sylph.tsv"
        run_cmd(f"sylph profile {ref_db} {self.filename} > {outfile}")

        for l in open(outfile):
            row = l.strip().split()
            print(row)
            