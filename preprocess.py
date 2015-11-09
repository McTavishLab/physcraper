import sys
from dendropy import DnaCharacterMatrix
infi=sys.argv[1]
outstub=sys.argv[2]
start=int(sys.argv[3])
stop=int(sys.argv[4])

orig = DnaCharacterMatrix.get(path=infi, schema="nexus")

prune = []
d = {}
for taxon, seq in orig.items():
	d[taxon.label] = seq.values()[start:stop]
  


dna = DnaCharacterMatrix.from_dict(d)

dna.write(path="{}.fas".format(outstub), schema="fasta")
