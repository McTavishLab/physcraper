from dendropy import Tree, DnaCharacterMatrix
import os
import configparser
import subprocess
import sys



runname="ascomycota"
config = configparser.ConfigParser()
config.read('/home/ejmctavish/projects/otapi/physcraper/config')

get_ncbi_taxonomy = config['ncbi.taxonomy']['get_ncbi_taxonomy']
ott_ncbi = config['ncbi.taxonomy']['ott_ncbi']
ncbi_dmp = config['ncbi.taxonomy']['ncbi_dmp']
MISSINGNESS_THRESH = 0.5

ncbi_to_ott = {}
fi =open(ott_ncbi)

#pickle meeeee
for lin in fi:
    lii= lin.split(",")
    ncbi_to_ott[int(lii[1])]=int(lii[0])

gi_ncbi_map = {}
if os.path.isfile("id_map.txt"):
    fi = open("id_map.txt")
    for lin in fi:
        gi_ncbi_map[int(lin.split(",")[0])]=lin.split(",")[1]


orig_seq = DnaCharacterMatrix.get(path="accs",schema="fasta")

#prune out identical sequences

mapped_taxon_ids=open("id_map.txt","a")
stops = []
for taxon, seq in orig_seq.items():
    gi = int(taxon.label.split('|')[1])
    if gi in gi_ncbi_map.keys():
        try:
            taxon.label = ncbi_to_ott[int(gi_ncbi_map[gi])]
        except:
            taxon.label = "ncbi_id_{}".format(int(gi_ncbi_map[gi]))
    else:
        try:
            ncbi_id = int(subprocess.check_output(["bash", get_ncbi_taxonomy, "{}".format(gi), "{}".format(ncbi_dmp)]).split('\t')[1])
            mapped_taxon_ids.write("{}, {}\n".format(gi, ncbi_id))
            gi_ncbi_map[gi] = ncbi_id
            try:
                taxon.label = ncbi_to_ott[int(gi_ncbi_map[gi])]
            except:
                taxon.label = "ncbi_id_{}".format(gi_ncbi_map[gi].strip())
        except:
            taxon.label = "gi_{}".format(gi)
            sys.stderr.write("no taxon id found for gi_{}".format(gi))
    stops.append(len(seq.values()))


stops.sort()
stop = stops[int(len(stops)/2)]

d = {}
for taxon, seq in orig_seq.items():
        d[str(taxon.label)] = seq.values()[:stop]
    

dna_orig = DnaCharacterMatrix.from_dict(d)

dna_taxa = [i for i in dna_orig.taxon_namespace]

tre_orig = Tree.get(path = "{}_random_resolve.tre".format("ascomycota"), schema = "newick",taxon_namespace=dna_orig.taxon_namespace)
""
treed_taxa = [i.taxon for i in tre_orig.leaf_nodes()]

tre_orig.prune_taxa(set(treed_taxa) - set(dna_taxa))

for taxon in set(dna_taxa) - set(treed_taxa):
	del d[taxon.label]

#####NEXT STEPS!!!

#make a function that doe sthis dumb shit in orig as well
dna_orig = DnaCharacterMatrix.from_dict(d)

tre_orig.write(path = "{}_orig_cut.tre".format(runname), schema = "newick", unquoted_underscores=True, suppress_edge_lengths=True)

dna_orig.write(path="{}_orig_cut.phy".format(runname), schema="phylip")
dna_orig.write(path="{}_orig_cut.fas".format(runname), schema="fasta")



p1 = subprocess.call(["papara", "-t","{}_orig_cut.tre".format(runname), "-s", "{}_orig_cut.phy".format(runname), "-n", "o_asc"]) 

