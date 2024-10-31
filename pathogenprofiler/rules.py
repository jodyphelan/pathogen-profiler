from typing import List, Union
from .models import Variant, Gene
import re
import logging

# Define all the necessary functions for the complete example

def parse_string(s: str) -> Union[bool, int, float, str]:
    """Parse string to most appropriate type"""
    if s.lower() == 'true':
        return True
    elif s.lower() == 'false':
        return False
    # if string consists of digits, convert to int
    elif s.isdigit():
        return int(s)
    # if string consists of digits and a decimal point, convert to float
    elif re.match(r'^\d+\.\d+$', s):
        return float(s)
    return s

def dsl_parse_genetic_object(query: str) -> dict:
    """
    Improved parsing of the DSL query to handle attribute keys and values correctly.

    Arguments
    ---------
    query : str
        The DSL query to parse.
    
    Returns
    -------
    dict
        The parsed query in a dictionary format.
    
    Examples
    --------
    >>> dsl_parse_genetic_object('Gene(name="gyrA",position="1234")')
    {'class': 'Gene', 'attributes': {'name': 'gyrA', 'position': 1234}}
    """
    pattern = r'(\w+)\((.*)\)'
    match = re.match(pattern, query)

    if not match:
        return None

    class_type, attributes = match.groups()
    attr_dict = {}

    for attr in attributes.split(','):
        key, value = attr.split('=')
        attr_dict[key.strip().strip('"').strip("'")] = parse_string(value.strip().strip('"').strip("'"))

    return {
        'class': class_type,
        'attributes': attr_dict
    }


def search_for_object(data_structure: List[Union[Variant,Gene]], query: str) -> List[dict]:
    """
    Search the data structure for a match based on the DSL query and return the matching object or None.

    Arguments
    ---------
    data_structure : List[Variant]
        The data structure to search.
    query : str
        The DSL query to search for.
    
    Returns
    -------
    dict|None
        The matching object or None if no match was found.
    
    Examples
    --------
    >>> from pathogenprofiler.data import Variant, Gene, generate_example_gene
    >>> genes = [generate_example_gene()]
    >>> search_for_object(genes, 'Gene(gene_id="Rv0667")')
    Gene(gene_id='Rv0667', type='functionally_normal', filter='pass', gene_name='rpoB', annotation=[{'type': 'drug_resistance', 'drug': 'rifampicin'}])
    >>> search_for_object(genes, 'Gene(gene_id="Rv0668")') == None
    True
    """
    parsed_query = dsl_parse_genetic_object(query)

    if not parsed_query:
        return None

    class_type = parsed_query['class']
    attributes = parsed_query['attributes']

    attribute_expansions = []
    if 'type' in attributes and attributes['type'] == 'lof':
        for t in ('frameshift_variant', 'stop_gained', 'transcript_ablation','feature_ablation'):
            tmp = attributes.copy()
            tmp['type'] = t
            attribute_expansions.append(tmp)
    else:
        attribute_expansions.append(attributes)

    matched_entries = []
    for attrib in attribute_expansions:
        for entry in data_structure:
            if isinstance(entry,Variant) and class_type=='Gene':
                continue
            if isinstance(entry,Gene) and class_type=='Variant':
                continue


            if all(vars(entry).get(key, None) == value for key, value in attrib.items()):
                matched_entries.append(entry)

    return matched_entries


def parse_flexible_dsl(rule: str) -> dict:
    """
    Parse the DSL query to handle actions between any two entities (genes or mutations) in an agnostic manner.

    Arguments
    ---------
    rule : str
        The DSL query to parse.

    Returns
    -------
    dict
        The parsed query in a dictionary format.
    
    Examples
    --------
    >>> parse_flexible_dsl('Variant(gene_name="gyrA") inactivates_resistance Variant(gene_name="gyrB")')
    {'source_query': 'Variant(gene_name="gyrA")', 'target_query': 'Variant(gene_name="gyrB")', 'action': 'inactivates_resistance'}
    """
    parts = rule.split()
    parts = parts[:3] + [' '.join(parts[3:])]
    if len(parts) < 3:
        return None
    if len(parts)==3:
        parts.append(None)
    source_query, action, target_query, note = parts

    return {
        'source_query': source_query,
        'target_query': target_query,
        'action': action,
        'note': note
    }


def execute_inactivates_resistance_flexible(data_structure: List[Variant], source_query: dict, target_query: dict, note: str, just_make_note=False) -> bool:
    """
    Execute the 'inactivates_resistance' action in a flexible manner where either genes or mutations
    can inactivate resistance on genes or mutations.
    """
    source_objects = search_for_object(data_structure, source_query)
    target_objects = search_for_object(data_structure, target_query)

    if source_objects and target_objects:
        for target_object in target_objects:
            for annotation in vars(target_object).get('annotation', []):
                if annotation.get('type') == 'drug_resistance':
                    if note:
                        annotation['note'] = note
                    if not just_make_note:
                        annotation['type'] = 'inactivated_drug_resistance'

        return True
    return False

def execute_make_interaction_note(data_structure: List[Variant], source_query: dict, target_query: dict) -> bool:
    """
    Execute the 'make_note' action in a flexible manner where either genes or mutations
    can inactivate resistance on genes or mutations.
    """
    source_objects = search_for_object(data_structure, source_query)
    target_objects = search_for_object(data_structure, target_query)
    if source_objects and target_objects:
        for target_object in target_objects:
            for annotation in vars(target_object).get('annotation', []):
                if annotation.get('type') == 'drug_resistance':
                    annotation['note'] = 'There is an interaction between %s and %s. Please check mutations for further details.' % (source_query, target_query)
        return True
    return False

def apply_rules(rules: List[str], genetic_objects: List[dict], just_make_note: bool=False) -> List[str]:
    """
    Apply the rules to the genetic objects.
    
    Arguments
    ---------
    rules : List[str]
        The rules to apply.
    genetic_objects : List[dict]
        The genetic objects to apply the rules to.
    
    Examples
    --------
    >>> from pathogenprofiler.data import Variant, Gene, Consequence, generate_example_gene
    >>> data = [
    ...     Variant(chrom='CU458896',pos=2345982,ref='T',alt='C',depth=100,freq=1,sv=False,filter='pass',gene_id='MAB_2297',gene_name='erm(41)',change='c.28T>C',nucleotide_change='c.28T>C',protein_change='p.Trp10Arg'),
    ...     Gene(gene_id="MAB_2297",type="functionally_normal",filter="pass",gene_name="erm(41)",annotation=[{"type":"drug_resistance","drug":"macrolides","literature":"10.1038/s41467-021-25484-9"}])
    ... ]
    >>> rules = ['Variant(gene_name="erm(41)", protein_change="p.Trp10Arg") inactivates_resistance Gene(gene_name="erm(41)")']
    >>> apply_rules(rules, data)
    >>> data[1].annotation[0]['type']=='inactivated_drug_resistance'
    True
    """
    actions = {
        'inactivates_resistance': execute_inactivates_resistance_flexible,
        'make_interaction_note': execute_make_interaction_note
    }
    rules_applied = []
    for rule in rules:
        logging.debug("Applying rule: %s" % rule)
        parsed_query = parse_flexible_dsl(rule)
        if parsed_query:
            action_executed = actions[parsed_query['action']](
                genetic_objects,
                parsed_query['source_query'], 
                parsed_query['target_query'],
                parsed_query['note'],
                just_make_note=just_make_note
            )
            if not action_executed:
                pass
            else:
                rules_applied.append(rule)


    return rules_applied
    