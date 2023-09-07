
import logging
from .utils import run_cmd
import csv
from uuid import uuid4

class SourmashSig:

    def __init__(self,filename,tmp_prefix=None):
        self.filename = filename
        self.tmp_prefix = tmp_prefix if tmp_prefix else str(uuid4())

    def filter(self,min_abundance=5,outfile=None):
        logging.info("Filtering sourmash sig")
        
        if outfile is None:
            outfile = self.filename.replace(".sig",".filtered.sig")
        
        run_cmd(f"sourmash sig filter -o {outfile} {self.filename} -m {min_abundance}")
        
        return SourmashSig(outfile,tmp_prefix=self.tmp_prefix)
    
    def search(self, ref_db, db_annotation, ani_threshold=0.95):
        logging.info("Searching sourmash sig")
        
        outfile = "%s" % self.tmp_prefix+".sourmash.csv"
        run_cmd(f"sourmash search -n 10 {self.filename} {ref_db} --containment --estimate-ani --ignore-abundance -o {outfile}")

        species = {}
        for row in csv.DictReader(open(db_annotation)):
            species[row["accession"]] = row["species"]

        results = []
        for row in csv.DictReader(open(outfile)):
            results.append({
                "accession": row["name"],
                "species":species[row["name"]],
                "ani":float(row["ani"])*100
            })
        results = [x for x in results if x["ani"]>=ani_threshold]
        return results[:10]

