# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 17:47:49 2018

@author: blubb

collection of functions to do little phylo tasks
    convert aln
    convert tree
"""




#### convert aln
input_fn = "/home/blubb/sync-TP-T470s/Senecioneae/its_sl2904.nex"
input_format="nexus"
output_fn = "/home/blubb/sync-TP-T470s/Senecioneae/its_sl.fasta"
output_format="fasta"


def convert_aln(input_fn, input_format, output_fn, output_format):
    """simple function to convert alignment files"""
    from Bio import SeqIO
    SeqIO.parse( input_fn, input_format)
    count = SeqIO.convert(input_fn, input_format, output_fn, output_format)
    print("Converted %i records" % count)

convert_aln(input_fn, input_format, output_fn, output_format)


####convert tree
tree_fn = "/home/blubb/sync-TP-T470s/Senecioneae/its_sl2904.nex.con.tre"
tree_format = "nexus"

output_tree_fn = "/home/blubb/sync-TP-T470s/Senecioneae/its_sl.nwk"
output_tree_format = "newick"


def convert_tree(tree_fn, tree_format, output_tree_fn,output_tree_format):
    """simple function to convert phylogenies to different formats"""
    import dendropy
    tree = dendropy.Tree.get(
        path=tree_fn,
        schema=tree_format)
    
    tree.write(
        path=output_tree_fn,
        schema=output_tree_format,
        )
    print("converted the tree")    
        
convert_tree(tree_fn, tree_format,output_tree_fn,tree_format)    