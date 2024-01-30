import pathogenprofiler as pp
import pytest
import logging

ref = "MTB-h37rv_asm19595v2-eg18.fa"
gff = "MTB-h37rv_asm19595v2-eg18.gff"


hgvs_mutations  = [
    {
        "Gene":"gid", 
        "Mutation":"c.615_628delACGACGTGGAAAGCinsGCGACGTGGAAAG",
        "target":"c.615_628delACGACGTGGAAAGCinsGCGACGTGGAAAG"
    },
    {
        "Gene":"gid", 
        "Mutation":"c.615_628delACGACGTGGAAAGC",
        "target":"c.615_628delACGACGTGGAAAGC"
    },
    {
        "Gene":"ethA", 
        "Mutation":"c.1057_1059del",
        "target":"c.1057_1059delACG"
    },
    {
        "Gene": "rpoB",
        "Mutation": "c.1282_1290del",
        "target": "c.1287_1295delGCTGAGCCA"
    },
    {
        "Gene":"ethA",
        "Mutation":"c.341del",
        "target":"c.341delA"
    },
    {
        "Gene": "embA" ,
        "Mutation": "c.-29_-28del",
        "target": "c.-29_-28delCT"
    },
    {
        "Gene": "alr",
        "Mutation": "c.-283_-280delCAAT",
        "target": "c.-283_-280delCAAT"
    },
    {
        "Gene": "ethA",
        "Mutation": "c.-1058_968del",
        "target": "c.-1058_968del"
    },
    {
        "Gene": "rpoB",
        "Mutation":"c.1296_1297insTTC",
        "target": "c.1297_1299dupTTC"
    },
    {
        "Gene": "pncA",
        "Mutation": "c.521_522insT",
        "target": "c.521_522insT"
    },
    {
        "Gene": "alr",
        "Mutation": "c.-283_-280delCAAT",
        "target": "c.-283_-280delCAAT"
    },
    {
        "Gene": "pncA",
        "Mutation": "c.7G>C",
        "target": "c.7G>C"
    },
    {
        "Gene": "rrs",
        "Mutation": "n.5G>A",
        "target": "n.5G>A"
    },
    {
        "Gene": "gid",
        "Mutation": "g.4408196C>G",
        "target": "c.7G>C",
    },
    {
        "Gene": "rrs",
        "Mutation": "n.-93C>T",
        "target": "n.-93C>T"
    },
    {
        "Gene": "rrs",
        "Mutation": "n.4_17dupTGTTTGGAGAGTTT",
        "target": "n.4_17dupTGTTTGGAGAGTTT"
    },
    {
        "Gene": "mmpR5",
        "Mutation": "c.-21_-20insTTC",
        "target": "c.-21_-20insTTC"
    },
    {
        "Gene": "gid",
        "Mutation": "c.114_115dupCC",
        "target": "c.114_115dupCC"
    },
    {
        "Gene": "gid",
        "Mutation": "c.248_250dupGGC",
        "target": "c.248_250dupGGC"
    },
    {
        "Gene": "dnaA",
        "Mutation": "c.-110_-109insTG",
        "target": "c.-110_-109insTG"
    }

]

todo = [
    {
        "Gene": "rrs",
        "Mutation": "n.-29_-28insATAC",
        "target": "n.-29_-28insATAC"
    },
]


converted_mutations = pp.get_snpeff_formated_mutation_list(hgvs_mutations,ref,gff,"Mycobacterium_tuberculosis_h37rv")




logging.debug(converted_mutations)

@pytest.mark.parametrize("mutation",hgvs_mutations)
def test_variant(mutation):
    gene = mutation["Gene"]
    change = mutation["Mutation"]
    target = mutation["target"]
    assert(target == converted_mutations[(gene,change)])

