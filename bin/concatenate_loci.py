#!/usr/bin/env python
import dendropy
import argparse
import os
import copy

from physcraper import aligntreetax

# Example using the data from https://github.com/McTavishLab/physcraper_example
#physcraper_run.py -s pg_238 -t tree109 -a physcraper_example/treebase_alns/pg_238tree109_18S.aln -db /branchinecta/shared/local_blast_db -nt 8 -as "nexus" -r -o pg_238_runs/18s
#physcraper_run.py -s pg_238 -t tree109 -a physcraper_example/treebase_alns/pg_238tree109_18S.aln -db /branchinecta/shared/local_blast_db -nt 8 -as "nexus" -r -o pg_238_runs/18s

#Examples
#python bin/concatenate_loci.py -d ../physcraper_runs/pg_238_runs/ -f concatenate -o test_concat_res



#python bin/concatenate_loci.py -d ../physcraper_runs/pg_238_runs/ -f astral -o test_concat_res
#java -jar astral.5.7.5.jar -i ../physcraper/test_concat_res/genetrees.new -a ../physcraper/test_concat_res/mappings.txt 


#ToDo - clean up dendropy part, write Charsets
# Generate taxon sets for SVD quratets inputs



parser = argparse.ArgumentParser()
parser.add_argument("-d","--locus_runs_folder", help="folder containing results directories from individual locus runs")
parser.add_argument("-o","--output", help="folder to write to")
parser.add_argument("-c","--concatenate_by", help="what aspect of otu info to concatenate by, eitehr OTTid or ...?")
parser.add_argument("-f","--format", help="output format", choices=["concatenate", "astral"])
parser.add_argument("-s","--schema", help="alignment file format schema", choices=["fasta", "nexus"], default="fasta")
parser.add_argument("-m","--include_missing", help="Where unevern numbers of sequences are available, concatenate with gaps", default = False)
## Include all seqs or just some???


args = parser.parse_args()



all_taxa = {}

#args.output = "test_concat_res"
if not os.path.exists(args.output):
    os.mkdir(args.output)

tax_labels = {}

otu_dict_dict = {}

runs = [f for f in os.listdir(args.locus_runs_folder) if os.path.isdir("{}/{}".format(args.locus_runs_folder, f))]
i = 0

all_loci = {}
for run in runs:
    print(run)
    i += 1
    fp = "{}/{}".format(args.locus_runs_folder, run)
    files = [f for f in os.listdir(fp)]
    if files:
        for file in files:
            if file.startswith('inputs_'):
                tag = file.split('.')[0].replace('inputs_', '')
        loc = "Locus{}".format(i)
        print(loc)
        locus = aligntreetax.generate_ATT_from_run(fp, start_files='output', tag=tag, run=False)
        all_loci[loc] = {}
        all_loci[loc]['data'] = locus
        all_loci[loc]['path'] = fp
        for seq in locus.aln:
            otu_id = seq.label
            taxon = locus.otu_dict[otu_id]['^ot:ottTaxonName']
            tax_labels[taxon] = locus.otu_dict[otu_id]['^ot:ottTaxonName']
            if taxon not in all_taxa:
                all_taxa[taxon] = {}
            if loc not in all_taxa[taxon]:
                all_taxa[taxon][loc] = []
            all_taxa[taxon][loc].append(otu_id)

nloci = i

#ds.read(path="pythonidae_cytb.fasta", schema="fasta", data_type="dna")
#ds.read(schema="pythonidae.mle.tre", "nexus", taxon_namespace=ds.taxon_namespaces[0])
#ds.write_to_path("pythonidae_combined.nex", "nexus")



if args.format == "concatenate":
    concat_dict = {}
    re_label_dict = {}
    for loc in all_loci:
        re_label_dict[loc] = {}

    taxa_in_all_loci = set(all_taxa.keys())

    for locus in all_loci:
        locusset = set([all_loci[locus]['data'].otu_dict[seq.label]['^ot:ottTaxonName'] for seq in all_loci[locus]['data'].aln])
        taxa_in_all_loci = taxa_in_all_loci.intersection(locusset)

    concat_names = []
    for taxon in taxa_in_all_loci:
        print(taxon)
        ntax_seq_max = max([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
        ntax_seq_min = min([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
        if args.include_missing:
            ntax_seq = ntax_seq_max
        else:
            ntax_seq = ntax_seq_min
        print(ntax_seq)
        for i in range(ntax_seq):
            concat_name = "{}_{}".format(tax_labels[taxon].replace(" ","_"), i)
            concat_dict[concat_name] = {}
            concat_names.append(concat_name)
            for loc in all_loci:
                print(loc)
                tip = all_taxa[taxon][loc][i]
                concat_dict[concat_name][loc] = tip
                re_label_dict[loc][tip] = concat_name
        #if args.include_missing == False:
        #    assert len(concat_dict[concat_name][loc]) == len(concat_names)


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
            else:
                seq.label = concat_name
        aln.remove_sequences(removal_list)
        ds.add(aln, taxon_namespace=ds.taxon_namespaces[0])


    ds.unify_taxon_namespaces()
    d_all = dendropy.DnaCharacterMatrix.concatenate(ds.char_matrices)
    d_all.write(path ="{}/concat.aln".format(args.output), schema = args.schema)


if args.format == "astral":
    genetrees = []
    taxset = {}
    for locus in all_loci:
        trecop = copy.deepcopy(all_loci[locus]['data'].tre)
        for tip in trecop.leaf_node_iter():
            tax = all_loci[locus]['data'].otu_dict[tip.taxon.label].get('^ot:ottTaxonName', "UNK")
            if tax not in taxset:
                taxset[tax] = []
            new_label = "{}_{}".format(locus, tip.taxon.label)
            tip.taxon.label = new_label
            taxset[tax].append(new_label)
        genetrees.append(trecop.as_string(schema = "newick"))
    
    fi = open("{}/genetrees.new".format(args.output), "w")
    for gtree in genetrees:
        gtree = gtree.replace('[&R] ', "")
        gtree = gtree.replace("'", "")
        fi.write(gtree)
    fi.close()
    
    mfi = open("{}/mappings.txt".format(args.output), "w")
    for tax in taxset:
        mfi.write(tax)
        mfi.write(":")
        mfi.write(",".join(taxset[tax]))
        mfi.write("\n")
    mfi.close()