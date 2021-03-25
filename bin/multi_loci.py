#!/usr/bin/env python
""" 
Example using the data from https://github.com/McTavishLab/physcraper_example
physcraper_run.py -s pg_238 -t tree109 -a physcraper_example/treebase_alns/pg_238tree109_18S.aln -db /branchinecta/shared/local_blast_db -nt 8 -as "nexus" -r -o pg_238_runs/18s
physcraper_run.py -s pg_238 -t tree109 -a physcraper_example/treebase_alns/pg_238tree109_18S.aln -db /branchinecta/shared/local_blast_db -nt 8 -as "nexus" -r -o pg_238_runs/18s

Examples:
python bin/multi_loci.py -d ../physcraper_runs/pg_238_runs/ -f concatenate -o test_concat_res


python bin/multi_loci.py -d ../physcraper_runs/pg_238_runs/ -f astral -o test_concat_res
java -jar astral.5.7.5.jar -i ../physcraper/test_concat_res/genetrees.new -a ../physcraper/test_concat_res/mappings.txt

"""
import argparse
import os
import copy
import dendropy

from physcraper import aligntreetax



parser = argparse.ArgumentParser()
parser.add_argument("-d", "--locus_runs_folder", help="folder containing results directories from individual locus runs")
parser.add_argument("-o", "--output", help="folder to write to")
#parser.add_argument("-c","--concatenate_by", help="what aspect of otu info to concatenate by, eitehr OTTid or ...?")
parser.add_argument("-f", "--format", help="output format", choices=["concatenate", "astral", "svdq"])
parser.add_argument("-s", "--schema", help="output alignment file format schema", choices=["fasta", "nexus"], default="fasta")
parser.add_argument("-m", "--include_missing", action="store_true", help="For concatenation, where uneven numbers of sequences are available, concatenate with gaps", default=False)
## Include all seqs or just some???

args = parser.parse_args()


assert args.output, "-o output folder is required"
assert args.locus_runs_folder, "-d input folder with completed runs is required"

if not os.path.exists(args.output):
    os.mkdir(args.output)


# Runs is the set of completed physcraper runs.
# They should all be in the folder given by -d aka locus_runs_folder

def setup_dicts(runs_path):
    """Inputs: a directory of completed physcraper runs
    Outputs: 
        all_loci: A dictionary of AlignTreeTaxon objects by locus
        all_taxa: a mapping of taxon names to tips in ATT objects
    """
    i = 0 #iterator to count runs (loci)
    all_loci = {} # dictionary of ATT objects
    all_taxa = {}
    runs = [f for f in os.listdir(runs_path) if os.path.isdir("{}/{}".format(runs_path, f))]
    for run in runs:
        print(run)
        i += 1
        fp = "{}/{}".format(runs_path, run)
        files = os.listdir(fp)
        if files:
            for file in files:
                if file.startswith('inputs_'):
                    tag = file.split('.')[0].replace('inputs_', '')
            loc = "Locus{}_{}".format(i, tag)
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
    return all_loci, all_taxa


def concatenate(all_loci, all_taxa, include_missing, gapchar="?"):
    """Concatenate loci from ATT objects stored in a dictionary
    Inputs:    
        loci: dictionary of ATT objects by locus
        taxa: dictionary of tips by taxon by locus
        include_missing: Boolean
            if False, drop taxa where a locus is missing
            if True, include taxa even if a locus is missing
        gapchar: Character to fill in for missing loci if include missing == T
                default = "?"
    Outputs:
        A dendropy data object
        """
    assert len(gapchar) == 1
    concat_dict = {}
    tax_map = {}
    for taxon in all_taxa:
        tax_map[taxon] = []
        ntax_seq_max = max([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
        # If we include missing data, this is how many represnetatives of this taxon we should end up with
        ntax_seq_min = min([len(all_taxa[taxon][locus]) for locus in all_taxa[taxon]])
        # If we do not include missing data, this is how many represnetatives of this taxon we should end up with
        if include_missing:
            ntax_seq = ntax_seq_max
        else:
            ntax_seq = ntax_seq_min
        #so after this, we have pruned down this set to only taxa that are present in all loci (stirnct misisng data condition)
        for i in range(ntax_seq):
            concat_name = "{}_{}".format(taxon.replace(" ", "_"), i)
            concat_dict[concat_name] = {}
            tax_map[taxon].append(concat_name)
#            concat_names.append(concat_name)
            for loc in all_loci:
                if i < len(all_taxa[taxon].get(loc, [])):
                    tip = all_taxa[taxon][loc][i]
                    concat_dict[concat_name][loc] = tip # tells me what sequence form this locus I will but in this concanated taxon
                else:
                    concat_dict[concat_name][loc] = "GAP"
    ds = dendropy.DataSet()
    taxon_namespace = dendropy.TaxonNamespace()
    ds.attach_taxon_namespace(taxon_namespace)
    re_label_dict = {locus:{} for locus in all_loci} # What are we calling each concatneated taxon, and what OTUs are part of it.
    no_missing_tax = []
    if not include_missing:
        for conc_tax in concat_dict:
            if 'GAP' not in concat_dict[conc_tax].values():
                no_missing_tax.append(conc_tax)
                for loc in concat_dict[conc_tax]:
                    tip = concat_dict[conc_tax][loc]
                    re_label_dict[loc][tip] = conc_tax
        for loc in all_loci:
            aln = copy.deepcopy(all_loci[loc]['data'].aln)
            removal_list = []
            for seq in aln:
                concat_name = re_label_dict[loc].get(seq.label, None)
                if concat_name is None:
                    removal_list.append(seq)
                else:
                    seq.label = concat_name
            aln.remove_sequences(removal_list)
            aln._label = loc
            ds.add(aln, taxon_namespace=ds.taxon_namespaces[0])

    if include_missing:
        for conc_tax in concat_dict:
            for loc in concat_dict[conc_tax]:
                tip = concat_dict[conc_tax][loc]
                re_label_dict[loc][tip] = conc_tax
        for loc in all_loci:
            aln = copy.deepcopy(all_loci[loc]['data'].aln)
            matched = set()
            for tax, seq in aln.items():
                seqlen = len(seq)
                concat_name = re_label_dict[loc].get(tax.label, None)
                matched.add(concat_name)
                tax.label = concat_name
            gap_seqs = {}
            for concat_name in concat_dict:
                if concat_name not in matched:
                    assert concat_dict[concat_name][loc] == "GAP"
                    gap_seqs[concat_name] = gapchar * seqlen
            gap_aln = dendropy.DnaCharacterMatrix.from_dict(gap_seqs, taxon_namespace=aln.taxon_namespace)
            aln.add_sequences(gap_aln)
            aln._label = loc
            ds.add(aln, taxon_namespace=ds.taxon_namespaces[0])
    ds.unify_taxon_namespaces()
    d_all = dendropy.DnaCharacterMatrix.concatenate(ds.char_matrices)
    return tax_map, concat_dict, d_all

def write_astral_inputs(all_loci, all_taxa, outputdir):
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
        genetrees.append(trecop.as_string(schema="newick"))

    fi = open("{}/genetrees.new".format(outputdir), "w")
    for gtree in genetrees:
        gtree = gtree.replace('[&R] ', "")
        gtree = gtree.replace("'", "")
        fi.write(gtree)
    fi.close()

    mfi = open("{}/mappings.txt".format(outputdir), "w")
    for tax in taxset:
        mfi.write(tax)
        mfi.write(":")
        mfi.write(",".join(taxset[tax]))
        mfi.write("\n")
    mfi.close()


def generate_taxpartition_str(tax_map):
    taxpart_str = "\ntaxpartition species =\n"
    for taxon in tax_map:
        taxpart_str = taxpart_str + "\t{}: ".format(taxon.replace(" ", "_").replace("-", "_")) +"'"+"' '".join(tax_map[taxon])+"',\n"
    taxpart_str = taxpart_str + ";\nend;"
    return taxpart_str

if args.format == "concatenate":
    loci, taxa = setup_dicts(args.locus_runs_folder)
    tax_map, concat_dict, d_all = concatenate(all_loci=loci,
                                              all_taxa=taxa,
                                              include_missing=args.include_missing,
                                              gapchar="?")
    d_all.write(path="{}/{}".format(args.output, 'concat.aln'), schema=args.schema)


if args.format == "svdq":
    loci, taxa = setup_dicts(args.locus_runs_folder)
    tax_map, concat_dict, d_all = concatenate(all_loci=loci,
                                              all_taxa=taxa,
                                              include_missing=args.include_missing,
                                              gapchar="?")
    nexstr = d_all.as_string(schema="nexus")
    taxpart_str = generate_taxpartition_str(tax_map)
    fi = open("{}/svdq.nex".format(args.output), "w")
    fi.write(nexstr[:-7]) #deletes "end;"
    fi.write(taxpart_str)
    fi.close()

if args.format == "astral":
    loci, taxa = setup_dicts(args.locus_runs_folder)
    write_astral_inputs(loci, taxa, args.output)

