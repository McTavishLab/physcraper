#!/usr/bin/env python
"""Concatenating across runs/genes"""

import pickle
import sys
import random
from dendropy import DnaCharacterMatrix



gene1 = pickle.load(open(sys.argv[1],"rb"))
gene2 = pickle.load(open(sys.argv[2],"rb"))

spp_to_otu1 = {}
spp_to_otu2 = {}


def get_spp_by_otu(physcraper_obj, spp_dict):
    avg_seqlen = sum(physcraper_obj.orig_seqlen)/len(physcraper_obj.orig_seqlen)
    seq_len_cutoff = avg_seqlen*physcraper_obj.seq_len_perc
    ok_otus = []
    for taxon, seq in physcraper_obj.aln.items():
        if len(seq.symbols_as_string().translate(None, "-?")) > seq_len_cutoff:
            ok_otus.append(taxon.label)
    seqlen = len(seq) #should all be same bc aligned
    for otu in physcraper_obj.otu_dict.keys():
        if otu in ok_otus:
            try:
                spp_name = physcraper_obj.otu_dict[otu]["^ot:ottId"]
                if spp_name:
                    spp_name = str(spp_name)
                    if spp_name not in spp_dict.keys():
                        spp_dict[spp_name] = []
                        spp_dict[spp_name].append(otu)
                    else:
                        spp_dict[spp_name].append(otu)
            except KeyError:
                pass

get_spp_by_otu(gene1, spp_to_otu1)
get_spp_by_otu(gene2, spp_to_otu2)

dellist = []

for spp in spp_to_otu1.keys():
    if spp not in spp_to_otu2:
        dellist.append(spp)


for spp in spp_to_otu2.keys():
    if spp not in spp_to_otu1:
        dellist.append(spp)

for item in dellist:
    try:
        del spp_to_otu1[item]
    except:
        del spp_to_otu2[item]

assert set(spp_to_otu1.keys()) == set(spp_to_otu2.keys())

def arbitrary_prune_fill(spp_dict, physcraper_obj):
    """prunes to 1 seq per spp, and fills in missing data for missing spp,
    in preparation for concanteneation, return dict to be made in char matrix"""
    aln_dict = {}
    tmp_dict = {}
    for taxon, seq in physcraper_obj.aln.items():
        aln_dict[taxon.label] = seq
    seqlen = len(seq) #should all be same bc aligned
    for spp_name in spp_dict.keys():
        try:
            otu = random.choice(spp_dict[spp_name])
            tmp_dict[spp_name] = aln_dict[otu]
        except KeyError:
            tmp_dict[spp_name] = "-" * seqlen
    return tmp_dict

aln1 = DnaCharacterMatrix.from_dict(arbitrary_prune_fill(spp_to_otu1, gene1))
aln2 = DnaCharacterMatrix.from_dict(arbitrary_prune_fill(spp_to_otu2, gene2), taxon_namespace = aln1.taxon_namespace)

concat = DnaCharacterMatrix.concatenate([aln1,aln2])
concat.write(path="concat.fas",
            schema="fasta")





#Open the two pyscraper objects
#Merge the alignements on OTT_ID?
#How to force/missing data ...


#Option 1: randomly select one seq from each ott ID.
#Option 2: Use all pairwise?
#Option 3: force mono phyly of spps?  grrrrrrrrrrrrrrrr