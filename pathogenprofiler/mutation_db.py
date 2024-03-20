from typing import List, Union, Set
import json
from .models import Variant, Gene, GenomePositionDepth
from copy import deepcopy
from collections import defaultdict
import re
from .hgvs import extract_numbers

class DictSet:
    """
    A class to create a set that works with dictionaries
    """
    def __init__(self,data: List[dict] = None):
        self.container = set()
        if data:
            for d in data:
                self.container.add(json.dumps(d))
    def add(self,data: Union[dict,list]) -> None:
        if isinstance(data,list):
            for d in data:
                self.container.add(json.dumps(d))
        elif isinstance(data,dict):
            self.container.add(json.dumps(data))

    def to_dict_list(self) -> List[dict]:
        return [json.loads(d) for d in self.container]



def db_compare(variants: List[Variant],db: dict) -> List[Union[Variant,Gene]]:
    annotated_results = deepcopy(variants)
    mutation_db = MutationDB(db)

    for var in annotated_results:
        mutation_db.annotate_variant(var)
        mutation_db.apply_lof_annotation(var)
    
    genes = mutation_db.get_functionally_normal_genes(annotated_results)
    annotated_results += genes

    return annotated_results

class MutationDB:
    """
    A class to store mutations in a data structure that can be queried
    Can be used to annotate mutations with data from a database
    """
    def __init__(self,db: dict):   
        self.db = {}
        self.genes = set()
        self.genome_pos2annotation = defaultdict(list)
        for gene in db:
            self.genes.add(gene)
            for var in db[gene]:
                self.db[(gene,var)] = db[gene][var]
                if self.db[(gene,var)]['genome_positions']:
                    for pos in self.db[(gene,var)]['genome_positions']:
                        for ann in self.db[(gene,var)]['annotations']:
                            # deep copy the annotation to avoid changing the original
                            ann = deepcopy(ann)
                            ann['gene'] = gene
                            ann['variant'] = var
                            self.genome_pos2annotation[(db[gene][var]['chromosome'],pos)].append(ann)
                        

    def annotate_variant(self,var: dict) -> List[dict]:
        """Annotate a variant with data from the database"""
        for csq in var.consequences:
            annotations = self.get_annotation(csq)
            csq.annotation = annotations

    def apply_lof_annotation(self,var):
        """Apply loss of function annotation to a variant"""
        for csq in var.consequences:
            for annotation in csq.annotation:
                if annotation['type']=='loss_of_function_variant':
                    csq.type = annotation['so_term']        
            

    def get_annotation(self,csq: dict) -> List[dict]:
        """Get the annotation for a consequence from the database"""
        db_var_match = DictSet()
        
        if csq.gene_id not in self.genes:
            return []
        
        key = (csq.gene_id,csq.nucleotide_change)
        if key in self.db:
            db_var_match.add(self.db[key]['annotations'])
        key = (csq.gene_id,csq.protein_change)
        if key in self.db:
            db_var_match.add(self.db[key]['annotations'])
        
        for t in csq.type.split("&"):
            key = (csq.gene_id,t)
            if key in self.db:
                db_var_match.add(self.db[(csq.gene_id,t)]['annotations'])
        if d:=self.check_for_so_wildcard(csq):
            db_var_match.add(d['annotations'])
        return  db_var_match.to_dict_list()
    
    def get_functionally_normal_genes(self, variants: List[dict]) -> List[str]:
        """Get all variants for a gene that are functionally normal"""
        genes_to_check =  [gene for gene,variant in self.db.keys() if variant=='functionally_normal']
        genes_to_return = []
        for gene in genes_to_check:
            intact = True
            for var in variants:
                for csq in var.consequences:
                    if csq.gene_id==gene:
                        if csq.type in ('loss_of_function_variant','stop_gained','frameshift_variant','feature_ablation','transcript_ablation'):
                            intact = False
            if intact:
                genes_to_return.append(gene)
        
        gene_objects = []
        for gene in genes_to_return:
            gene_objects.append(Gene(
                gene_id=gene,
                annotation=self.db[(gene,'functionally_normal')]['annotations'],
                type="functionally_normal"
            ))
        
        return gene_objects
    

    def get_gene_variants(self,gene: str) -> List[str]:
        """Get all variants for a gene"""
        return [v for g,v in self.db if g==gene]
    
    def check_for_so_wildcard(self,csq: dict):
        """Check if the variant is in the database with a wildcard SO term"""
        for var in self.get_gene_variants(csq.gene_id):
            for t in csq.type.split('&'):
                r = re.search(f"{t}_([pcn])\.(\d+)_(\d+)",var)
                if r: 
                    context = r.group(1)
                    positions = set(range(int(r.group(2)),int(r.group(3))+1))
                    if context=="p":
                        change = csq.protein_change
                    elif context=="c":
                        change = csq.nucleotide_change
                    elif context=="n":
                        change = csq.nucleotide_change
                    affected_positions = extract_affected_positions(change)
                    if len(positions.intersection(affected_positions))>0:
                        return self.db[(csq.gene_id,var)]
                    
    def annotate_missing_positions(self, positions: List[GenomePositionDepth]) -> None:
        """
        Annotate missing positions with data from the database
        
        Arguments
        ---------
        positions: List[PositionDepth]
            A list of positions to annotate
        
        Returns
        -------
        List[MissingPos]
            A list of annotated positions
        """
        annotated_positions = []
        for pos in positions:
            if (pos.chrom,pos.pos) in self.genome_pos2annotation:
                pos.annotation = self.genome_pos2annotation[(pos.chrom,pos.pos)]
        return annotated_positions

def extract_affected_positions(change: str) -> Set[int]:
    """Extract the affected positions from a variant annotation with a range of positions"""
    numbers = extract_numbers(change)
    return set(range(numbers[0],numbers[-1]+1))
    