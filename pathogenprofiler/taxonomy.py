
import logging
from typing import List
from .utils import run_cmd
import csv
from uuid import uuid4
import os
from .models import TaxonomicHit, Species


def combine_species_abundance(matches: List[dict]) -> List[dict]:
    """
    Combine species abundance from multiple hits
    Parameters
    ----------
    matches : List[dict]
        List of taxonomic hits
    
    Returns
    -------
    List[dict]
        List of combined species hits with relative abundance
    """
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
    """Base class for sketching and taxonomic classification"""
    def __init__(self,filename, tmp_prefix=None):
        self.filename = filename
        self.tmp_prefix = tmp_prefix if tmp_prefix else str(uuid4())

    def get_species_hits(self, ref_db: str, db_annotation: str, *args, **kwargs) -> List[Species]:
        """
        Get species hits with annotations
        Parameters
        ----------
        ref_db : str
            Path to the reference database
        db_annotation : str
            Path to the database annotation CSV file
        *args, **kwargs
            Additional arguments for the classify method
            See subclass implementations for details
        Returns
        -------
        List[Species]
            List of species hits with annotations
        """
        sequence_hits = self.classify(ref_db, *args, **kwargs)

        accession_data = {}
        for row in csv.DictReader(open(db_annotation)):
            accession_data[row["accession"]] = row

        species_hits = []
        for hit in sequence_hits:
            data = hit.model_dump()
            data.update(accession_data.get(hit.accession,{}))
            species_hits.append(Species(**data))

        combined_species_hits = combine_species_abundance(species_hits)
        return combined_species_hits
    
    def classify(self, ref_db: str, *args, **kwargs) -> List[TaxonomicHit]:
        raise NotImplementedError("Subclasses must implement classify method")
    


class SourmashSig(Sketch):

    def __init__(self,filename: str,tmp_prefix: str = None):
        super().__init__(filename, tmp_prefix=tmp_prefix)

    def classify(self, ref_db: str, intersect_bp: int=500000,f_match_threshold: float=0.1) -> List[TaxonomicHit]:
        """
        Classify using sourmash
        
        Parameters
        ----------
        ref_db : str
            Path to the sourmash reference database
        intersect_bp : int
            Minimum intersect base pairs to consider a hit
        f_match_threshold : float
            Minimum f_match to consider a hit

        Returns
        -------
        List[TaxonomicHit]
            List of taxonomic hits
        """
        outfile = "%s" % self.tmp_prefix+".sourmash.csv"
        run_cmd(f"sourmash gather {self.filename} {ref_db} -o {outfile}")

    
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
    
    def __repr__(self):
        return "SourmashSig(%s)" % self.filename


class SylphSketch(Sketch):
    
    def __init__(self,filename,tmp_prefix=None):
        super().__init__(filename, tmp_prefix=tmp_prefix)

    def classify(self, ref_db: str, threads: int = 1) -> List[TaxonomicHit]:
        """
        Classify using sylph
        
        Parameters
        ----------
        ref_db : str
            Path to the sylph reference database directory
        threads : int
            Number of threads to use

        Returns
        -------
        List[TaxonomicHit]
            List of taxonomic hits
        """
        outfile = "%s" % self.tmp_prefix+".sylph.tsv"
        run_cmd(f"sylph profile -t {threads} {ref_db}/* {self.filename} > {outfile}")

        hits = []
        for row in csv.DictReader(open(outfile), delimiter="\t"):
            hit = TaxonomicHit(
                prediction_method="sylph",
                accession=row["Genome_file"].split("/")[-1].replace(".fasta","").replace(".fna","").replace(".fa",""),
                ani=round(float(row["Adjusted_ANI"]),2),
                abundance=float(row["Taxonomic_abundance"])
            )
            hits.append(hit)
        return hits
    
    def __repr__(self):
        return "SylphSketch(%s)" % self.filename