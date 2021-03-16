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
#parser.add_argument("-c","--concatenate_by", help="what aspect of otu info to concatenate by, eitehr OTTid or ...?")
parser.add_argument("-f","--format", help="output format", choices=["concatenate", "astral"])
parser.add_argument("-s","--schema", help="alignment file format schema", choices=["fasta", "nexus"], default="fasta")
parser.add_argument("-m","--include_missing", help="Where uneven numbers of sequences are available, concatenate with gaps", default = False)
## Include all seqs or just some???

args = parser.parse_args()

if not os.path.exists(args.output):
    os.mkdir(args.output)


# Track what taxa are present across all loci
all_taxa = {}

# find otu_ids by taxon
tax_labels = {}

# Runs is the set of completed physcraper runs. 
# They should all be in the folder given by -d aka locus_runs_folder
runs = [f for f in os.listdir(args.locus_runs_folder) if os.path.isdir("{}/{}".format(args.locus_runs_folder, f))]


i = 0 #iterator to count runs (loci)
all_loci = {} # dictionary of ATT objects
for run in runs:
    print(run)
    i += 1
    fp = "{}/{}".format(args.locus_runs_folder, run)
    files = [f for f in os.listdir(fp)]
    if files:
        for file in files:
            if file.startswith('inputs_'):
                tag = file.split('.')[0].replace('inputs_', '')
        loc = "Locus{}_{}".format(i, tag)
        print(loc)
        locus = aligntreetax.generate_ATT_from_run(fp, start_files='output', tag=tag, run=False) #build ATT opject for each run
        all_loci[loc] = {}
        all_loci[loc]['data'] = locus #all_loci[loc]['data'] has an alignemnt, and otuu_dictionary, and a tree
        all_loci[loc]['path'] = fp
        for seq in locus.aln:
            otu_id = seq.label
            taxon = locus.otu_dict[otu_id]['^ot:ottTaxonName']
            if taxon not in all_taxa:
                all_taxa[taxon] = {}
            if loc not in all_taxa[taxon]:
                all_taxa[taxon][loc] = []
            all_taxa[taxon][loc].append(otu_id)

nloci = i


if args.format == "concatenate":
    concat_dict = {}
    re_label_dict = {} # What are we calling each concatneated taxon, and what OTUs are part of it.
    for loc in all_loci:
        re_label_dict[loc] = {}

  
    for taxon in all_taxa:
        ntax_seq_max = max([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]]) 
        # If we include missing data, this is how many represnetatives of this taxon we should end up with
        ntax_seq_min = min([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
        # If we do not include missing data, this is how many represnetatives of this taxon we should end up with
        if args.include_missing:
            ntax_seq = ntax_seq_max
        else:
            ntax_seq = ntax_seq_min
        #so after this, we have pruned down this set to only taxa that are present in all loci (stirnct misisng data condition)
        for i in range(ntax_seq):
            concat_name = "{}_{}".format(taxon.replace(" ","_"), i)
            concat_dict[concat_name] = {}
#            concat_names.append(concat_name)
            for loc in all_loci:
                if i <= len(all_taxa[taxon][loc]):
                    tip = all_taxa[taxon][loc][i]
                    concat_dict[concat_name][loc] = tip # tells me what sequence form this locus I will but in this concanated taxon
                    re_label_dict[loc][tip] = concat_name
                else:
                    concat_dict[concat_name][loc] = "GAP"
        #if args.include_missing == False:
        #    assert len(concat_dict[concat_name][loc]) == len(concat_names)

    ds = dendropy.DataSet()
    taxon_namespace = dendropy.TaxonNamespace()
    ds.attach_taxon_namespace(taxon_namespace)

    for loc in all_loci:
        aln = copy.deepcopy(all_loci[loc]['data'].aln)
        if args.include_missing == False:
            removal_list = []
            for seq in aln:
                concat_name = re_label_dict[loc].get(seq.label, None)
                    if concat_name == None:
                        removal_list.append(seq)
                    else:
                        seq.label = concat_name
            aln.remove_sequences(removal_list)
        if args.include_missing == True:
            matched = set()
            for seq in aln:
                seqlen = len(seq.)
                concat_name = re_label_dict[loc].get(seq.label, None)
                assert concat_name
                matched.add(concat_name)
            for concat_name in concat_dict:
                if concat_name not in matched:
                    assert concat_dict[concat_name][loc] = "GAP"
                newseq = "-" * seqlen
                
                #Add newseq to 



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