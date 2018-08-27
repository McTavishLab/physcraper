# # altered from https://github.com/zyxue/ncbitax2lin/blob/master/ncbitax2lin.py

import re
import os
import gzip
import multiprocessing
# import argparse
import logging
import sys
import pandas as pd
# from utils import timeit, backup_file

nodes = None
names = None


def strip(str_):
    '''    :param str_: a string
    '''
    return str_.strip()

class parser:
    """read in databases from ncbi to connect species name with the taxonomic identifier 
    and the corresponding hirachical information. It provides a much faster way to get those information.
    We use those files to get independent from web requests to find those information 
    (the implementation of it in BioPython was not really reliable).
    Nodes includes the hierachical infortmation, names the scientific names and ID's.
    The files need to be updated regularly, best way to always do it when a new blast database was loaded
    """
    def __init__(self, names_file, nodes_file):
        # self.nodes = self.load_nodes(nodes_file)
        # self.names = self.load_names(names_file)
        self.names_file = names_file
        self.nodes_file = nodes_file
        self.initialize()

    def initialize(self):
        """ The data itself are not stored in __init__, as then the information will be pickled (which results in gigantic pickle file sizes).
        Instead everytime the function is loaded after loading a pickle file, it will be 'initialized'.
        """
        global nodes
        nodes = self.load_nodes(self.nodes_file)
        global names
        names = self.load_names(self.names_file)

    def load_nodes(self, nodes_file):
        """ load nodes.dmp and convert it into a pandas.DataFrame
        """ 
        # print(nodes_file)
        assert os.path.exists(nodes_file)
        df = pd.read_csv(nodes_file, sep='|', header=None, index_col=False,
                         names=[
                             'tax_id',
                             'parent_tax_id',
                             'rank',
                             'embl_code',
                             'division_id',
                             'inherited_div_flag',
                             'genetic_code_id',
                             'inherited_GC__flag',
                             'mitochondrial_genetic_code_id',
                             'inherited_MGC_flag',
                             'GenBank_hidden_flag',
                             'hidden_subtree_root_flag',
                             'comments'
                         ])

        # To get rid of flanking tab characters
        df['rank'] = df['rank'].apply(strip)
        df['embl_code'] = df['embl_code'].apply(strip)
        df['comments'] = df['comments'].apply(strip)
        return df

    def load_names(self, names_file):
        """ Load names.dmp and convert it into a pandas.DataFrame
        """ 
        assert os.path.exists(names_file)

        df = pd.read_csv(names_file, sep='|', header=None, index_col=False,
                         names=[
                             'tax_id',
                             'name_txt',
                             'unique_name',
                             'name_class'
                         ])
        df['name_txt'] = df['name_txt'].apply(strip)
        df['unique_name'] = df['unique_name'].apply(strip)
        df['name_class'] = df['name_class'].apply(strip)

        sci_df = df[df['name_class'] == 'scientific name']
        sci_df.reset_index(drop=True, inplace=True)
        return sci_df


    def get_downtorank_id(self, tax_id, downtorank="species"):
        """ Recursive function to find the parent id of a taxon as defined by downtorank.
        """ 
        # print(tax_id)
        # print(type(tax_id))
        # print(downtorank)
        if nodes is None:
            self.initialize()

        if type(tax_id) != int:
            sys.stdout.write("WARNING: tax_id {} is no integer. Will convert value to int\n".format(tax_id))
            tax_id = int(tax_id)
        # assert type(tax_id) is int()


        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == downtorank:
            return tax_id
        else:
            parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            return self.get_downtorank_id(parent_id, downtorank)

    def get_name_from_id(self, tax_id):
        """ Find the scientific name for a given ID.
        """ 
        if names is None:
            self.initialize() 
        return names[names["tax_id"] == tax_id]["name_txt"].values[0].replace(" ", "_")

