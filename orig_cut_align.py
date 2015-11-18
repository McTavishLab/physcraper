from dendropy import Tree, DnaCharacterMatrix
import os
import configparser
import subprocess

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
        	taxon.label = "ncbi_id_{}".format(gi_ncbi_map[gi])
    else:
        ncbi_id = int(subprocess.check_output(["bash", get_ncbi_taxonomy, "{}".format(gi), "{}".format(ncbi_dmp)]).split('\t')[1])
        mapped_taxon_ids.write("{}, {}\n".format(gi, ncbi_id))
        gi_ncbi_map[gi] = ncbi_id
        try:
        	taxon.label = ncbi_to_ott[int(gi_ncbi_map[gi])]
        except:
        	taxon.label = "ncbi_id_{}".format(gi_ncbi_map[gi])
    stops.append(len(seq.values()))


stop = sum(stops)/len(stops)

d = {}
for taxon, seq in orig_seq.items():
        d[taxon.label] = seq.values()[start:stop]
    

dna_cut = DnaCharacterMatrix.from_dict(d)

for taxon in d.taxon_namespace:
        if taxon.label.replace("_"," ") in map_dict:
            if  map_dict[taxon.label.replace("_"," ")] == {}:
                taxon.label = taxon.label.replace("/","_") # it's legal nexus, but RAxML chokes. Probably other chars this is true of as well...
            else:
                if map_dict[taxon.label.replace("_"," ")]['^ot:ottId'] not in mapped_ids: #Can't have two tips with same name. Need better alternative tho!
                    mapped_ids.add(map_dict[taxon.label.replace("_"," ")]['^ot:ottId'])
                    taxon.label = str(map_dict[taxon.label.replace("_"," ")]['^ot:ottId'])
                else:
                    taxon.label = taxon.label.replace("/","_")
        else:
            sys.sterr.write("taxon label problem")


  
dna_cut = DnaCharacterMatrix.from_dict(d)
tre.prune_taxa_with_labels(exclude)

tre.write(path = "{}_cut.tre".format(runname), schema = "newick", unquoted_underscores=True, suppress_edge_lengths=True)

dna_cut.write(path="{}_aln_ott_cut.phy".format(runname), schema="phylip")
dna_cut.write(path="{}_aln_ott_cut.fas".format(runname), schema="fasta")

