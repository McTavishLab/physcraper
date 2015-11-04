
from dendropy import DnaCharacterMatrix
infi=sys.argv[1]
outstub=sys.argv[2]
start=int(sys.argv[3])
stop=int(sys.argv[4])

orig = DnaCharacterMatrix.get(path="M4058.nexorg", schema="nexus")

prune = []
d = {}
for taxon, seq in orig.items():
	d[taxon.label] = seq.values()[start:stop]
	if len(seq.values()[start:stop].translate(None, "-?")) == 0:
            prune.append(taxon.label)
  
for tax in prune:
	del d[tax]


dna = DnaCharacterMatrix.from_dict(d)

tre.prune_taxa_with_labels(prune)

dna.write(path="{}.fas".format(outstub), schema="fasta")


 
  
    