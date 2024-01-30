from .models import Variant, Consequence

def set_change(var: Variant) -> None:
    """
    Set the change field for a variant
    
    Parameters
    ----------
    var : Variant
        The variant to set the change field for
    """
    protein_csqs = ["missense_variant","stop_gained"]
    var.change = var.protein_change if var.type in protein_csqs else var.nucleotide_change
    return var

def select_most_relevant_csq(csqs):
    rank = ["transcript_ablation","exon_loss_variant","frameshift_variant","large_deletion","start_lost","disruptive_inframe_deletion","disruptive_inframe_insertion","stop_gained","stop_lost","conservative_inframe_deletion","conservative_inframe_insertion","initiator_codon_variant","missense_variant","non_coding_transcript_exon_variant","upstream_gene_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_variant","stop_retained_variant","splice_region_variant","synonymous_variant"]

    ranked_csq = sorted(csqs,key=lambda x: min([rank.index(y) if y in rank else 999 for y in x['type'].split("&")]))
    return ranked_csq[0]


def select_csq(var: Variant) -> None:
    """Select the most relevant consequence for a variant"""
    annotated_csq = []
    for csq in var.consequences:
        if len(csq.annotation)>0:
            annotated_csq.append(csq)
    if len(annotated_csq)==0:
        csq = select_most_relevant_csq(var.consequences)
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
        alternate_consequences = [json.dumps(x) for x in d["consequences"]]
        alternate_consequences.remove(json.dumps(csq))
        alternate_consequences = [json.loads(x) for x in alternate_consequences]
    del d["consequences"]
    d.update(csq)
    d["alternate_consequences"] = alternate_consequences
    d = set_change(d)

    