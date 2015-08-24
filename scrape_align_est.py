#!/usr/bin/env python
import sys
from dendropy import Tree, DnaCharacterMatrix
from peyotl import gen_otu_dict, iter_node
from peyotl.manip import iter_trees
from peyotl.phylesystem.phylesystem_umbrella import Phylesystem
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree,  PhyloSchema
from peyotl.api import APIWrapper
api_wrapper = APIWrapper()
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO, Entrez
import re
import subprocess
import time

'''
#LSU ASC tree example
study_id = "pg_873"
tree_id = "tree1679"
seqaln = "mat1.fas"
mattype="fasta"
runname="asc_test"
'''

study_id=sys.argv[1]
tree_id=sys.argv[2]
seqaln=sys.argv[3]
mattype=sys.argv[4]
runname=sys.argv[5]
firstrun=sys.argv[6]

remote=1

#Fixed values
E_VALUE_THRESH = 0.04
ott_ncbi="/home/ejmctavish/projects/otapi/physcraper/ott_ncbi" #TODO config file
Entrez.email = "ejmctavish@gmail.com"





ott_to_ncbi = {}
ncbi_to_ott = {}
fi =open(ott_ncbi)
for lin in fi:
    lii= lin.split(",")
    ott_to_ncbi[int(lii[0])]=int(lii[1])
    ncbi_to_ott[int(lii[1])]=int(lii[0])


fi.close()

ottids = []

if firstrun:
    sys.stdout.write("First run through")
    phy = Phylesystem()
    n = phy.return_study(study_id)[0]
    api_wrapper.study.get(study_id,tree=tree_id)
    ##This is a weird way to get the ingroup node, but I need the OTT ids anyhow.

    #Now grab the same tree with the orginal lablels
    newick = extract_tree(n, tree_id, PhyloSchema('newick', output_nexml2json = '1.2.1', content="tree", tip_label="ot:originalLabel"))
    newick = newick.replace(" ", "_") #UGH
    d = DnaCharacterMatrix.get(path=seqaln,
                               schema=mattype)
    d.taxon_namespace.is_mutable = True
    tre = Tree.get(data=newick,
                    schema="newick",
                    taxon_namespace=d.taxon_namespace)

    # get all of the taxa associated with tips of the tree, and make sure that
    #   they include all of the members of the data's taxon_namespace...
    treed_taxa = [i.taxon for i in tre.leaf_nodes()]
    if len(treed_taxa) != len(d.taxon_namespace):
        missing = [i.label for i in d.taxon_namespace if i not in treed_taxa]
        emf = 'Some of the taxa in the alignment are not in the tree. Missing "{}"\n'
        em = emf.format('", "'.join(missing))
        raise ValueError(em)

    #Get original label to ott id mapping to match alignement to tree 
    ps = PhyloSchema('nexson', content='otumap', version='1.2.1')
    map_dict = map_dict=ps.convert(n)

    for taxon in d.taxon_namespace:
        if taxon.label.replace("_"," ") in map_dict:
            if  map_dict[taxon.label.replace("_"," ")] == {}:
                pass
            else:
                taxon.label = str(map_dict[taxon.label.replace("_"," ")]['^ot:ottId'])
                ottids.append(taxon.label)
        else:
            pass

    tre.resolve_polytomies()
    tre.write(path = "{}_random_resolve.tre".format(runname), schema = "newick", unquoted_underscores=True, suppress_edge_lengths=True)

    d.write(path="{}_aln_ott.phy".format(runname), schema="phylip")
    d.write(path="{}_aln_ott.fas".format(runname), schema="fasta")
else:
    d = DnaCharacterMatrix.get(path="{}_aln_ott.fas".format(runname),
                               schema="fasta",
                               preserve_underscores=True)
    d.taxon_namespace.is_mutable = False
    tre = Tree.get(path="{}_stream.tre".format(runname),
                    schema="newick",
                    taxon_namespace=d.taxon_namespace)
    d.write(path="{}_aln_ott.phy".format(runname), schema="phylip")


    # get all of the taxa associated with tips of the tree, and make sure that
    #   they include all of the members of the data's taxon_namespace.
    # This shouldn't fail as we repeat through the loop!

    for taxon in d.taxon_namespace:
        try:
            if int(taxon.label) in ott_to_ncbi.keys():
                ottids.append(taxon.label)
            else:
                pass
        except:
            pass


assert(len(ottids) > 2)


#This section grabs the MRCA node and blasts for seqs that are desc from that node
mrca_node = tree_of_life.mrca(ott_ids=ottids, wrap_response=True)
print("mrca_node found")


if remote:
    equery = "txid{}[orgn]".format(ott_to_ncbi[mrca_node.nearest_taxon.ott_id])
    for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
        record.seq._data = record.seq._data.replace("-","") # blast gets upset about too many gaps from aligned file
        print("blasting seq {}".format(i))
        result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"),  entrez_query=equery)
        save_file = open("{}_{}.xml".format(runname,i), "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()

else:
    for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
        output_handle = open("query.fasta", "w")
        SeqIO.write(record, output_handle, "fasta")
        output_handle.close()
        blastx_cline =  NcbiblastxCommandline(cmd='blastn', out="{}_{}.xml".format(runname,i), outfmt=5, query="query.fasta", db='nt', evalue=E_VALUE_THRESH)
        stdout, stderr = blastx_cline()


gi_to_ncbi = {}
new_seqs={}

for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
    result_handle = open("{}_{}.xml".format(runname,i))
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                   new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                   gi_to_ncbi[int(alignment.title.split('|')[1])] = ''

##SET MINIMUM SEQUENCE LENGTH

for gi in gi_to_ncbi.keys():
    gi_to_ncbi[gi] = int(subprocess.check_output(["bash", "get_ncbi_taxonomy.sh",  "{}".format(gi)]).split('\t')[1])


newseqs = "{}.fasta".format(runname)
fi = open(newseqs,'w')
print("writing out sequences")
for gi in gi_to_ncbi:
        try:
            ott_id = ncbi_to_ott[gi_to_ncbi[gi]]
        except:
            print("ncbi taxon ID {} not in ott".format(gi_to_ncbi[gi]))
            continue
        if ott_id not in ottids: # only adding seqs we don't have taxa for
                ottids.append(ott_id)
                fi.write(">{}\n".format(ncbi_to_ott[gi_to_ncbi[gi]]))
                fi.write("{}\n".format(new_seqs[gi]))
                print("success {}".format(ott_id))
        if ott_id in ottids: 
                print("ncbi taxon ID {} already in tree".format(gi_to_ncbi[gi]))


fi.close()


p1 = subprocess.call(["papara", "-t","{}_random_resolve.tre".format(runname), "-s", "{}_aln_ott.phy".format(runname), "-q",  newseqs, "-n", "extended"]) 
                          #run RAXML EPA on the alignments

#placement

p2 = subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-f", "v", "-s", "papara_alignment.extended", "-t","{}_random_resolve.tre".format(runname), "-n", "{}_PLACE".format(runname)])


placetre = Tree.get(path="RAxML_labelledTree.{}_PLACE".format(runname),
                schema="newick",
                preserve_underscores=True)
placetre.resolve_polytomies()
placetre.write(path = "{}_place_resolve.tre".format(runname), schema = "newick", unquoted_underscores=True)


#Full run with starting tree from placements

subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-s", "papara_alignment.extended", "-t","{}_place_resolve.tre".format(runname), "-p", "1", "-n", "{}".format(runname)])


e = DnaCharacterMatrix.get(path="papara_alignment.extended",
                           schema="phylip")

e.taxon_namespace.is_mutable = False # Names should be the same at this point.

newtre = Tree.get(path="RAxML_bestTree.{}".format(runname),
                schema="newick",
                taxon_namespace=e.taxon_namespace,
                preserve_underscores =True)


newtre.write(path = "{}_stream.tre".format(runname), schema = "newick", unquoted_underscores=True)
d.write(path="{}_aln_ott.fas".format(runname), schema="fasta")