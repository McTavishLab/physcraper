#!/usr/bin/env python
import dendropy
import argparse
import os
import copy


# Example using the data from https://github.com/McTavishLab/physcraper_example
#physcraper_run.py -s pg_238 -t tree109 -a physcraper_example/treebase_alns/pg_238tree109_18S.aln -db /branchinecta/shared/local_blast_db -nt 8 -as "nexus" -r -o pg_238_runs/18s
#physcraper_run.py -s pg_238 -t tree109 -a physcraper_example/treebase_alns/pg_238tree109_18S.aln -db /branchinecta/shared/local_blast_db -nt 8 -as "nexus" -r -o pg_238_runs/18s






from physcraper import aligntreetax
parser = argparse.ArgumentParser()
parser.add_argument("-d","--loci_runs_folder", help="folder containing results directories from individual locus runs")
parser.add_argument("-o","--output", help="folder to write to")
parser.add_argument("-c","--concatenate_by", help="what aspect of otu info to concatenate by, eitehr OTTid or ...?")
parser.add_argument("-m","--include_missing", help="Where unevern numbers of sequences are available, concatenate with gaps", default = False)
## Include all seqs or just some???


args = parser.parse_args()


ds = dendropy.DataSet()

# Set it up to manage all data under a single taxon namespace.
# HIGHLY RECOMMENDED!
taxon_namespace = dendropy.TaxonNamespace()
ds.attach_taxon_namespace(taxon_namespace)


all_taxa = {}

args.output = "test_concat_res"
if not os.path.exists(args.output):
    os.mkdir(args.output)

tax_labels = {}

otu_dict_dict = {}

runs = [f for f in os.listdir(args.loci_runs_folder)]
i = 0

all_loci = {}
for run in runs:
    print(run)
    i += 1
    fp = "{}/{}".format(args.loci_runs_folder, run)
    files = [f for f in os.listdir(fp)]
    for file in files:
        if file.startswith('inputs_'):
            tag = file.split('.')[0].replace('inputs_', '')
    loc = "Locus{}".format(i)
    print(loc)
    locus = aligntreetax.generate_ATT_from_run(fp, start_files='output', tag=tag, run=False)
    all_loci[loc] = {}
    all_loci[loc]['data'] = locus
    all_loci[loc]['path'] = fp
    for tip in locus.otu_dict:
        taxon = locus.otu_dict[tip]['^ot:ottId']
        tax_labels[taxon] = locus.otu_dict[tip]['^ot:ottTaxonName']
        if taxon not in all_taxa:
            all_taxa[taxon] = {}
        if loc not in all_taxa[taxon]:
            all_taxa[taxon][loc] = []
        all_taxa[taxon][loc].append(tip)

nloci = i


#
#ds.read(path="pythonidae_cytb.fasta", schema="fasta", data_type="dna")
#ds.read(schema="pythonidae.mle.tre", "nexus", taxon_namespace=ds.taxon_namespaces[0])
#ds.write_to_path("pythonidae_combined.nex", "nexus")


concat_dict = {}
re_label_dict = {}
for loc in all_loci:
    re_label_dict[loc] = {}


for taxon in all_taxa:
    ntax_seq_max = max([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
    ntax_seq_min = min([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
    if ntax_seq_min < len(all_loci):
        ntax_seq_min = 0
    if args.include_missing:
        ntax_seq = ntax_seq_max
    else:
        ntax_seq = ntax_seq_min
    for i in range(ntax_seq):
        concat_name = "{}_{}".format(tax_labels[taxon], i)
        concat_dict[concat_name] = {}
        for loc in all_loci:
            print(loc)
            if i < len(all_taxa[taxon].get(loc,[])):
                tip = all_taxa[taxon][loc][i]
                concat_dict[concat_name][loc] = tip
                re_label_dict[loc][tip] = concat_name

            else:
                concat_dict[concat_name][loc] = None


ds = dendropy.DataSet()
taxon_namespace = dendropy.TaxonNamespace()
ds.attach_taxon_namespace(taxon_namespace)


for loc in all_loci:
    aln = copy.deepcopy(all_loci[loc]['data'].aln)
    removal_list = []
    for seq in aln:
        concat_name = re_label_dict[loc].get(seq.label, None)
        if concat_name == None:
            removal_list.append(seq)
            print("remove")
            print(seq)
        else:
            seq.label = concat_name
    aln.remove_sequences(removal_list)
    ds.add(aln, taxon_namespace=ds.taxon_namespaces[0])


ds.unify_taxon_namespaces()

d_all = dendropy.DnaCharacterMatrix.concatenate(ds.char_matrices)
