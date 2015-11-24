from dendropy import Tree, DnaCharacterMatrix
import sys


d = {}

query_seq = DnaCharacterMatrix.get(path="ascomycota.fasta",schema="fasta")

def seq_dict_build(seq, label, seq_dict):
    new_seq = seq.symbols_as_string().replace("-","")
    for tax in seq_dict.keys():
        inc_seq = seq_dict[tax].symbols_as_string().replace("-","")
        if len(inc_seq) > len(new_seq):
            if inc_seq.find(new_seq) != -1:
                sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax))
                return
        else:
            if new_seq.find(inc_seq) != -1:
                del d[tax]
                d[label] = seq
                sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax, label, tax))
                return
    print (".")
    d[label] = seq
    return


for taxon, seq in query_seq.items():
    if len(seq.values()) > 800:
        seq_dict_build(seq, taxon.label, d)
    else:
        sys.stdout.write("*")


cull = DnaCharacterMatrix.from_dict(d)
cull.write(path="query_cull.fas", schema="fasta")