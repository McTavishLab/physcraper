"""
Function to plot an output tree and save it as pdf
"""

import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo

def plot(tree_file, pdf_file, fig_size = (7, 7), font_size = 12):
    # read the tree_file
    input_tree = Phylo.read(tree_file, "newick")
    # plot the tree
    # define the font size, it can be a variable later:
    matplotlib.rc('font', size= font_size)
    # set the size of the figure
    fig = plt.figure(figsize = fig_size, dpi=300)
    # alternatively
    # fig.set_size_inches(10, 40)
    # set width of the plot, maybe? explore this later
    # How else can we set this
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(input_tree, axes=axes, do_show = False)
    # save as pdf
    plt.savefig(pdf_file, dpi=300)
