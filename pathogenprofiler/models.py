from typing import List, Union, Optional, Generator
from pydantic import BaseModel, Field
from statistics import median

def generate_example_consequence():
    """Generate an example consequence"""
    return Consequence(
        gene_id='Rv0667',
        gene_name='rpoB',
        feature_id='CCP43410',
        type='missense_variant',
        nucleotide_change='c.1349C>T',
        protein_change='p.Ser450Leu',
        annotation=[{'type':'drug_resistance','drug':'rifampicin'}]
    )

def generate_example_variant(forward_reads=50,reverse_reads=50):
    """Generate an example variant"""
    return Variant(
        chrom='Chromosome',
        pos=761155,
        ref='C',
        alt='T',
        depth=forward_reads+reverse_reads,
        freq=0.5,
        forward_reads=forward_reads,
        reverse_reads=reverse_reads,
        sv=False,
        sv_len=False,
        consequences=[
            Consequence(
                gene_id='Rv0667',
                gene_name='rpoB',
                feature_id='CCP43410',
                type='missense_variant',
                nucleotide_change='c.1349C>T',
                protein_change='p.Ser450Leu',
                annotation=[{'type':'drug_resistance','drug':'rifampicin'}]
            )
        ]
    )

def generate_example_dr_variant():
    """Generate an example drug resistant variant"""
    return DrVariant(
        chrom='Chromosome',
        pos=761155,
        ref='C',
        alt='T',
        depth=100,
        freq=0.5,
        forward_reads=50,
        reverse_reads=50,
        sv=False,
        sv_len=0,
        consequences=[
            Consequence(
                gene_id='Rv0667',
                gene_name='rpoB',
                feature_id='CCP43410',
                type='missense_variant',
                nucleotide_change='c.1349C>T',
                protein_change='p.Ser450Leu',
                annotation=[{'type':'drug_resistance','drug':'rifampicin'}]
            )
        ],
        drugs=[{'drug':'rifampicin'}]
    )

def generate_example_gene():
    """Generate an example gene"""
    return Gene(
        gene_id='Rv0667',
        type='functionally_normal',
        gene_name='rpoB',
        annotation=[{'type':'drug_resistance','drug':'rifampicin'}]
    )

def generate_example_dr_gene():
    """Generate an example drug resistant gene"""
    return DrGene(
        gene_name='rpoB',
        gene_id='Rv0667',
        annotation=[],
        drugs=[{'drug':'rifampicin'}]
    )

class Consequence(BaseModel):
    """
    A consequence of a variant
    
    Attributes
    ----------
    gene_id : str
        The gene id
    gene_name : str
        The gene name
    feature_id : str
        The feature id
    type : str
        The type of variant
    nucleotide_change : str
        The nucleotide change
    protein_change : str
        The protein change
    annotation : List[dict]
        A list of annotations

    Examples
    --------
    >>> from pathogenprofiler import Consequence, generate_example_consequence
    >>> csq = generate_example_consequence()
    >>> csq
    Consequence(gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
    """
    gene_id: str
    gene_name: Union[str,None]
    feature_id: Union[str,None]
    type: str
    nucleotide_change: str
    protein_change: Union[str,None]
    sequence_hgvs: Union[str,None] = None
    annotation: List[dict] = Field(default_factory=list)

    def causes_drug_resistance(self, drug: str = None) -> bool:
        """
        Check if variant causes drug resistance

        Parameters
        ----------
        drug : str
            The drug to check for

        Returns
        -------
        bool
            True if variant causes drug resistance, False otherwise
        
        Examples
        --------
        >>> from pathogenprofiler import Consequence, generate_example_consequence
        >>> csq = generate_example_consequence()
        >>> csq.causes_drug_resistance()
        True
        >>> csq.causes_drug_resistance('rifampicin')
        True
        >>> csq.causes_drug_resistance('isoniazid')
        False
        """
        for ann in self.annotation:
            if ann['type']=='drug_resistance':
                if drug is None:
                    return True
                if ann['drug']==drug:
                    return True
        return False


class Variant(BaseModel):
    """
    A variant

    Attributes
    ----------
    chrom : str
        The chromosome
    pos : int
        The position
    ref : str
        The reference allele
    alt : str
        The alternate allele
    depth : int
        The depth
    freq : float
        The frequency
    sv : bool
        True if variant is a structural variant, False otherwise
    filter : str
        The filter string (pass, soft_fail or hard_fail)
    forward_reads : int
        The number of forward reads
    reverse_reads : int
        The number of reverse reads
    sv_len : int
        The length of the structural variant
    consequences : List[Consequence]
        A list of consequences
    gene_id : str
        The gene id
    gene_name : str
        The gene name
    feature_id : str
        The feature id
    type : str
        The type of variant
    change : str
        The change
    nucleotide_change : str
        The nucleotide change
    protein_change : str
        The protein change
    annotation : List[dict]
        A list of annotations

    Examples
    --------
    >>> from pathogenprofiler import Variant, generate_example_variant
    >>> var = generate_example_variant()
    >>> var
    Variant(chrom='Chromosome', pos=761155, ref='C', alt='T', depth=100, freq=0.5, forward_reads=50, reverse_reads=50, sv=False, sv_len=0, consequences=[Consequence(gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])], gene_id=None, gene_name=None, feature_id=None, type=None, change=None, nucleotide_change=None, protein_change=None, annotation=[])
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    depth: int
    freq: float
    sv: bool
    filter: Optional[str] = None
    forward_reads: Optional[int] = None
    reverse_reads: Optional[int] = None
    sv_len: Optional[int] = None
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    feature_id: Optional[str] = None
    type: Optional[str] = None
    change: Optional[str] = None
    nucleotide_change: Optional[str] = None
    protein_change: Optional[str] = None
    annotation: List[dict] = Field(default_factory=list)
    consequences: List[Consequence] = Field(default_factory=list)

    def select_most_relevant_csq(self):
        """
        Select the most relevant consequence for a variant
        
        Examples
        --------
        >>> from pathogenprofiler import Variant, generate_example_variant
        >>> var = generate_example_variant()
        >>> var.select_most_relevant_csq()
        Consequence(gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
        """

        rank = ["transcript_ablation","exon_loss_variant","frameshift_variant","large_deletion","start_lost","disruptive_inframe_deletion","disruptive_inframe_insertion","stop_gained","stop_lost","conservative_inframe_deletion","conservative_inframe_insertion","initiator_codon_variant","missense_variant","non_coding_transcript_exon_variant","upstream_gene_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_variant","stop_retained_variant","splice_region_variant","synonymous_variant"]

        ranked_csq = sorted(self.consequences,key=lambda x: min([rank.index(y) if y in rank else 999 for y in x.type.split("&")]))
        return ranked_csq[0]

    def set_default_csq(self):
        """
        Set the default consequence for a variant
        
        Examples
        --------
        >>> from pathogenprofiler import Variant, generate_example_variant
        >>> var = generate_example_variant()
        >>> var.set_default_csq()
        >>> var
        Variant(chrom='Chromosome', pos=761155, ref='C', alt='T', depth=100, freq=0.5, forward_reads=50, reverse_reads=50, sv=False, sv_len=0, consequences=[Consequence(gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])], gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', change='p.Ser450Leu', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
        """
        annotated_csq = []
        for csq in self.consequences:
            if len(csq.annotation)>0:
                annotated_csq.append(csq)
        if len(annotated_csq)==0:
            rank = ["transcript_ablation","exon_loss_variant","frameshift_variant","large_deletion","start_lost","disruptive_inframe_deletion","disruptive_inframe_insertion","stop_gained","stop_lost","conservative_inframe_deletion","conservative_inframe_insertion","initiator_codon_variant","missense_variant","non_coding_transcript_exon_variant","upstream_gene_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_variant","stop_retained_variant","splice_region_variant","synonymous_variant"]
            ranked_csq = sorted(self.consequences,key=lambda x: min([rank.index(y) if y in rank else 999 for y in x.type.split("&")]))
            csq = ranked_csq[0]
        elif len(annotated_csq)==1:
            csq = annotated_csq[0]
        else:
            chosen_annotation = None
            for csq in annotated_csq:
                if csq.causes_drug_resistance()==True:
                    chosen_annotation = csq
                    break
            if chosen_annotation:
                csq = chosen_annotation
            else:
                csq = annotated_csq[0]

        vars(self).update(vars(csq))
        protein_csqs = ["loss_of_function_variant","missense_variant","stop_gained","start_lost"]
        self.change = self.protein_change if self.type in protein_csqs else self.nucleotide_change    
    
    def set_gene_name(self, gene_names: dict) -> None:
        """
        Reset gene names for a variant
        
        Parameters
        ----------
        gene_names : dict
            A dictionary of gene names
        
        Examples
        --------
        >>> from pathogenprofiler import Variant, generate_example_variant
        >>> var = generate_example_variant()
        >>> var.set_default_csq()
        >>> var.set_gene_name({'Rv0667':'rpoB'})
        >>> var
        Variant(chrom='Chromosome', pos=761155, ref='C', alt='T', depth=100, freq=0.5, forward_reads=50, reverse_reads=50, sv=False, sv_len=0, consequences=[Consequence(gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])], gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', change='p.Ser450Leu', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
        """
        self.gene_name = gene_names[self.gene_id]
        for csq in self.consequences:
            csq.gene_name = gene_names.get(csq.gene_id,csq.gene_name)

    def convert_to_dr_element(self) -> None:
        """
        Convert a variant to a drug resistant variant
        
        Examples
        --------
        >>> from pathogenprofiler import Variant, generate_example_variant
        >>> var = generate_example_variant()
        >>> var.set_default_csq()
        >>> var.convert_to_dr_element()
        >>> var
        DrVariant(chrom='Chromosome', pos=761155, ref='C', alt='T', depth=100, freq=0.5, forward_reads=50, reverse_reads=50, sv=False, sv_len=0, consequences=[Consequence(gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])], gene_id='Rv0667', gene_name='rpoB', feature_id='CCP43410', type='missense_variant', change='p.Ser450Leu', nucleotide_change='c.1349C>T', protein_change='p.Ser450Leu', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}], drugs=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
        """
        if any([csq.causes_drug_resistance() for csq in self.consequences]):
            self.__class__ = DrVariant
            self.drugs = []
            for ann in self.annotation:
                if ann['type']=='drug_resistance':
                    self.drugs.append(ann)

    def get_str(self):
        """
        Return string representation of variant

        Examples
        --------
        >>> from pathogenprofiler import DrVariant, generate_example_dr_variant
        >>> var = generate_example_dr_variant()
        >>> var.get_str()
        'Chromosome 761155 C>T (0.50)'
        """
        return "%s %s (%.2f)" % (self.gene_name,self.change,self.freq)

    def get_annotation_value(self, annotation_type: str, key: str) -> List[str]:
        """
        Get annotation value
        
        Parameters
        ----------
        annotation_type : str
            The annotation type
        key : str
            The key to return
        
        Returns
        -------
        str
            The annotation value
        
        Examples
        --------
        >>> from pathogenprofiler import DrVariant, generate_example_dr_variant
        >>> var = generate_example_dr_variant()
        >>> var.get_annotation_value('drug_resistance','drug')
        'rifampicin'
        """
        return [ann[key] for ann in self.annotation if ann['type']==annotation_type]

    def __lt__(self, other) -> bool:
        return self.get_str() < other.get_str()
    
    def __hash__(self) -> int:
        return self.model_dump_json().__hash__()

class DrVariant(Variant):
    drugs: List[dict] = Field(default_factory=list)
    
    def get_drugs(self):
        """
        Return list of drugs associated with variant or gene
        
        Examples
        --------
        >>> from pathogenprofiler import DrVariant, generate_example_dr_variant
        >>> var = generate_example_dr_variant()
        >>> var.get_drugs()
        ['rifampicin']
        """
        return [x['drug'] for x in self.drugs]


class Gene(BaseModel):
    """
    A gene
    
    Attributes
    ----------
    gene_id : str
        The gene id
    type: str
        The sequence ontology term for the variant (usually functionally_normal)
    gene_name : str
        The gene name
    annotation : List[dict]
        A list of annotations

    Examples
    --------
    >>> from pathogenprofiler import Gene, generate_example_gene
    >>> gene = generate_example_gene()
    >>> gene
    Gene(gene_id='Rv0667', type='functionally_normal', gene_name='rpoB', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
    """
    gene_id: str
    type: str = 'functionally_normal'
    filter: str = 'pass'
    gene_name: Optional[str] = None
    annotation: List[dict] = Field(default_factory=list)

    def set_gene_name(self, gene_names: dict) -> None:
        """
        Reset gene names for a gene
        
        Parameters
        ----------
        gene_names : dict
            A dictionary of gene names
        
        Examples
        --------
        >>> from pathogenprofiler import Gene, generate_example_gene
        >>> gene = generate_example_gene()
        >>> gene.set_gene_name({'Rv0667':'rpoB'})
        >>> gene
        Gene(gene_id='Rv0667', type='functionally_normal', gene_name='rpoB', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
        """
        self.gene_name = gene_names[self.gene_id]
    
    def convert_to_dr_element(self) -> None:
        """
        Convert a gene to a drug resistant gene
        
        Examples
        --------
        >>> from pathogenprofiler import Gene, generate_example_gene
        >>> gene = generate_example_gene()
        >>> gene.convert_to_dr_element()
        >>> gene
        DrGene(gene_id='Rv0667', type='functionally_normal', gene_name='rpoB', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}], drugs=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
        """
        if any([ann['type']=='drug_resistance' for ann in self.annotation]):
            self.__class__ = DrGene
            self.drugs = []
            for ann in self.annotation:
                if ann['type']=='drug_resistance':
                    self.drugs.append(ann)
    def get_str(self):
        """
        Return string representation of variant

        Examples
        --------
        >>> from pathogenprofiler import Gene, generate_example_gene
        >>> gene = generate_example_gene()
        >>> gene.get_str()
        'rpoB (functionally_normal)'
        """
        return "%s (%s)" % (self.gene_name,self.type)

    def __lt__(self, other) -> bool:
        return True
    
class DrGene(Gene):
    """
    A gene that causes drug resistance
    
    Attributes
    ----------
    gene_id : str
        The gene id
    type: str
        The sequence ontology term for the variant (usually functionally_normal)
    gene_name : str
        The gene name
    annotation : List[dict]
        A list of annotations
    drugs : List[dict]
        A list of drugs associated with the gene

    Examples
    --------
    >>> from pathogenprofiler import DrGene
    >>> gene = DrGene(
    ...    gene_name='rpoB',
    ...    gene_id='Rv0667',
    ...    type='functionally_normal',
    ...    annotation=[{'type':'drug_resistance','drug':'rifampicin'}],
    ...    drugs=[{'drug':'rifampicin'}]
    ... )
    >>> gene
    DrGene(gene_id='Rv0667', type='functionally_normal', gene_name='rpoB', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}], drugs=[{'drug': 'rifampicin'}])
    """
    drugs: List[dict] = Field(default_factory=list)
    
    def get_drugs(self):
        """
        Return list of drugs associated with variant or gene
        
        Examples
        --------
        >>> from pathogenprofiler import DrGene
        >>> gene = DrGene(
        ...    gene_name='rpoB',
        ...    gene_id='Rv0667',
        ...    annotation=[],
        ...    type='functionally_normal',
        ...    drugs=[{'drug':'rifampicin'}]
        ... )
        >>> gene.get_drugs()
        ['rifampicin']
        """
        return [x['drug'] for x in self.drugs]
    def get_str(self):
        """
        Return string representation of variant

        Examples
        --------
        >>> from pathogenprofiler import DrGene
        >>> gene = DrGene(
        ...    gene_name='rpoB',
        ...    gene_id='Rv0667',
        ...    annotation=[],
        ...    type='functionally_normal',
        ...    drugs=[{'drug':'rifampicin'}]
        ... )
        >>> gene.get_str()
        'rpoB (functionally_normal)'
        """
        return "%s (resistance_gene:%s)" % (self.gene_name,self.type)


class GenomePosition(BaseModel):
    chrom: str
    pos: int

    def __hash__(self) -> int:
        return (self.chrom, self.pos).__hash__()
    
    def __lt__(self, other) -> bool:
        if self.chrom != other.chrom:
            raise ValueError(f"Cannot compare positions on different chromosomes: {self.chrom} and {other.chrom}")
        return self.pos < other.pos

class BarcodePosition(BaseModel):
    id: str
    chrom: str
    pos: int
    target_allele_count: int
    other_allele_count: int
    all_allele_count: int
    target_allele_percent: float


class GenomePositionDepth(BaseModel):
    """
    Class storing information about a missing position
    
    Attributes
    ----------
    chromosome : str
        The chromosome
    position : int
        The position
    depth : int
        The depth
    annotation : List[dict]
        A list of annotations

    Examples
    --------
    >>> from pathogenprofiler import GenomePosition
    >>> missing_pos = GenomePosition(
    ...    chromosome='Chromosome',
    ...    position=761155,
    ...    depth=1,
    ...    annotation=[{'type':'drug_resistance','drug':'rifampicin'}]
    ... )
    >>> missing_pos
    GenomePosition(chromosome='Chromosome', position=761155, depth=1, annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
    """
    chrom: str
    pos: int
    depth: int = None
    annotation: List[dict] = []

class GenomeRange(BaseModel):
    chrom: str
    start: int
    end: int

    def iter_positions(self) -> Generator[GenomePosition,None,None]:
        for p in range(self.start,self.end):
            yield GenomePosition(chrom=self.chrom,pos=p)

    def __contains__(self,pos: GenomePosition):
        if pos.chrom==self.chrom and pos.pos >= self.start and pos.pos <= self.end:
            return True
        else:
            return False
        
    def __hash__(self) -> int:
        return (self.chrom, self.start, self.end).__hash__()

class TargetQC(BaseModel):
    """
    Class storing information about the coverage of a target

    Attributes
    ----------
    target : str
        The target region
    percent_depth_pass : float
        The percent depth pass
    median_depth : float
        The median depth

    Examples
    --------
    >>> from pathogenprofiler import TargetQC
    >>> target_qc = TargetQC(
    ...    target='rpoB',
    ...    percent_depth_pass=100,
    ...    median_depth=100
    ... )
    >>> target_qc
    TargetQC(target='rpoB', percent_depth_pass=100.0, median_depth=100)
    """

    target: str
    percent_depth_pass: float
    median_depth: float

class SequenceQC(BaseModel):
    def get_target_median_depth(self) -> Optional[float]:
        if hasattr(self,'target_qc'):
            return median([x.median_depth for x in self.target_qc])
        else:
            return None
        
    def get_reads_mapped(self) -> Optional[int]:
        if hasattr(self,'num_reads_mapped'):
            return self.num_reads_mapped
        else:
            return None
        
    def get_percent_reads_mapped(self) -> Optional[float]:
        if hasattr(self,'percent_reads_mapped'):
            return self.percent_reads_mapped
        else:
            return None

class FastqQC(SequenceQC):
    num_sequences: int
    num_bases: int


class FastaQC(SequenceQC):
    """
    Class storing information about the quality of a fasta file

    Attributes
    ----------
    num_sequences : int
        The number of sequences
    num_bases : int
        The number of bases
    n50 : int
        The n50
    target_qc : List[TargetQC]
        A list of target qc objects

    Examples
    --------
    >>> from pathogenprofiler import FastaQC, TargetQC
    >>> fasta_qc = FastaQC(
    ...    num_sequences=20,
    ...    num_bases=4000000,
    ...    n50=200000,
    ...    target_qc=[
    ...        TargetQC(
    ...            target='rpoB',
    ...            percent_depth_pass=100,
    ...            median_depth=1
    ...        )
    ...    ]
    ... )
    >>> fasta_qc   
    FastaQC(num_sequences=20, num_bases=4000000, n50=200000, target_qc=[TargetQC(target='rpoB', percent_depth_pass=100.0, median_depth=1)])
    """
    num_sequences: int
    num_bases: int
    n50: int
    target_qc: List[TargetQC]


class BamQC(SequenceQC):
    """
    Class storing information about quality of the data in a bam file

    Attributes
    ----------
    percent_reads_mapped : float
        The percent reads mapped
    num_reads_mapped : int
        The number of reads mapped
    target_median_depth : float
        The median depth of targets
    genome_median_depth : int
        The median depth of the genome
    target_qc : List[TargetQC]
        A list of target qc objects

    Examples
    --------
    >>> from pathogenprofiler import BamQC, TargetQC, GenomePosition
    >>> bam_qc = BamQC(
    ...    percent_reads_mapped=100,
    ...    num_reads_mapped=100,
    ...    target_median_depth=100,
    ...    genome_median_depth=100,
    ...    target_qc=[
    ...        TargetQC(
    ...            target='rpoB',
    ...            percent_depth_pass=100,
    ...            median_depth=100
    ...        )
    ...    ]
    ... )
    >>> bam_qc
    BamQC(percent_reads_mapped=100.0, num_reads_mapped=100, target_median_depth=100, genome_median_depth=100, target_qc=[TargetQC(target='rpoB', percent_depth_pass=100.0, median_depth=100)])
    """

    percent_reads_mapped: float
    num_reads_mapped: int
    target_median_depth: float
    genome_median_depth: Optional[float]
    target_qc: List[TargetQC]
    missing_positions: List[GenomePositionDepth] = []


class VcfQC(SequenceQC):
    """
    Class storing information about quality of the data in a vcf file

    Attributes
    ----------
    num_variants : int
        The number of variants
    """
    total_variants: int
    


class SourmashSpeciesInfo(BaseModel):
    species: str
    accession: str
    ani: float



    
class BarcodeResult(BaseModel):
    """
    Class storing information about a prediction using a mutation barcode

    Attributes
    ----------
    id : str
        The id
    frequency: float
        The frequency
    info : list
        A list of info

    """

    id: str
    frequency: float
    info: list = []
    support: List[BarcodePosition]







class Species(BaseModel):
    """
    Class storing information about a species
    
    Attributes
    ----------
    species : str
        The species
    prediction_info : dict
        The prediction info
        
    Examples
    --------
    >>> from pathogenprofiler import Species
    >>> species = Species(
    ...    species='Mycobacterium abscessus',
    ...    prediction_info=None
    ... )
    >>> species
    Species(species='Mycobacterium abscessus', prediction_info=None)
    """
    species: str
    ani: Optional[float] = None
    abundance: Optional[float] = None
    relative_abundance: Optional[float] = None
    coverage: Optional[float] = None
    accession: Optional[str] = None
    ncbi_organism_name: Optional[str] = None
    prediction_method: Optional[str] = None
    notes: Optional[List[str]] = []


class TaxonomicHit(BaseModel):
    prediction_method: str = None
    accession: str
    abundance: Optional[float] = None
    ani: Optional[float] = None
    num_reads: Optional[int] = None

class SpeciesPrediction(BaseModel):
    """
    Class storing information about a species
    
    Attributes
    ----------
    prediction_method : str
        The prediction method
    species : List[Species]
        A list of species objects
    species_db : dict
        The species database

    """
    taxa: List[Species] = []
    qc_fail_taxa: List[Species] = []
    species_db: dict = {}

    def get_species_str(self):
        """
        Return the species
        
        Examples
        --------
        >>> from pathogenprofiler import SpeciesPrediction, Species
        >>> species_prediction = SpeciesPrediction(
        ...    prediction_method='user_specified',
        ...    species=[
        ...        Species(
        ...            species='Mycobacterium abscessus',
        ...            prediction_info=None
        ...        )
        ...    ]
        ... )
        >>> species_prediction.get_species()
        Mycobacterium abscessus
        """
        return ';'.join([x.species for x in self.species])




class ProfileResult(BaseModel):
    id: str
    software_version: str
    input_data_source: str 
    database: dict
    pipeline_software: list = Field(default_factory=list)
    variants: List[Union[Variant,Gene]]

class BamProfileResult(ProfileResult):
    qc: BamQC

class VcfProfileResult(ProfileResult):
    input_data_source: str = 'vcf'