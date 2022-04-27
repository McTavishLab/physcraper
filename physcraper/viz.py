"""
Function to plot a phylogenetic tree and save it as pdf
"""

import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo
import pylab
import csv


def plot_tree(tree_file,
              tree_file_format = "newick",
              pdf_file = "tree.pdf",
              pdf_size = (7, 7),
              font_size = 12,
              node_label_color = "blue",
              axis = "off",
              otu_info_file = "None"):
    # read the tree_file
    input_tree = Phylo.read(tree_file, tree_file_format)
    # plot the tree
    # define the font size, it can be a variable later:
    matplotlib.rc('font', size = font_size, weight = "bold")
    # define node labels color
    matplotlib.rc('text', color = node_label_color)
    # set the size of the figure
    fig = plt.figure(figsize = pdf_size, dpi = 300)
    # alternatively
    # fig.set_size_inches(10, 40)
    # set width of the plot, maybe? explore this later
    # How else can we set this
    axes = fig.add_subplot(1, 1, 1)
    # define the type of axis, default ot "off"
    pylab.axis(axis)
    # set tip label colors
    # read otu info file:
    # if otu_info_file != "None":
    #     otu_info =
    # input_tree = Phylo.read(t_file, "newick")
    # for tip in input_tree.get_terminals():
    #     print(tip.name)
    tip_label_colors = {"Brachychiton_acerifolius_otu376453": "r"}
    # draw the tree
    Phylo.draw(input_tree,
               axes=axes,
               do_show = False,
               label_colors = tip_label_colors)
    # save as pdf
    plt.savefig(pdf_file, dpi=300)
