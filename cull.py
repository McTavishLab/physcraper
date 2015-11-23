from dendropy import Tree, DnaCharacterMatrix
import sys


d = {}

query_seq = DnaCharacterMatrix.get(path="ascomycota.fasta",schema="fasta")

seqs = []

def contains(small, big):
    for i in xrange(1 + len(big) - len(small)):
        if small == big[i:i+len(small)]:
            return 1
    return 0


def seq_dict_build(seq, label, seq_dict):
    for tax in seq_dict.keys():
        inc_seq = seq_dict[tax]
        if len(inc_seq) > len(seq):
            if contains(seq, ex_seq):
                sys.stdout.write("seq {} is subsequence of {}, not added".format(label, tax))
                return
        else:
            if contains(ex_seq, seq):
                del d[tax]
                d[label] = seq
                sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed".format(label, tax, label, tax))
                return
    print (".")
    d[label] = seq
    return


for taxon, seq_in in query_seq.items():
    seq = seq_in.values()
    if len(seq) > 800:
        seq_dict_build(seq, taxon.label, d)
    else:
        sys.stdout.write("*")


cull = DnaCharacterMatrix.from_dict(d)
cull.write(path="query_cull.fas".format(runname), schema="fasta")