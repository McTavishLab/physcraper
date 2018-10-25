"""uses ncbi databases to easily retrieve taxonomic information.

parts are altered from https://github.com/zyxue/ncbitax2lin/blob/master/ncbitax2lin.py
"""

import os
import sys
import pandas as pd


_DEBUG_MK = 0

def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)

debug("Current ncbi_parser version number: 10172018.0")

nodes = None
names = None


def strip(str_):
    """ Strips of blank characters from string in pd dataframe.
    """
    return str_.strip()


def load_nodes(nodes_file):
    """ Loads nodes.dmp and converts it into a pandas.DataFrame.
    Contains the information about the taxonomic hierarchy of names.
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


def load_names(names_file):
    """ Loads names.dmp and converts it into a pandas.DataFrame.
    Includes only names which are accepted as scientific name by ncbi.
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


def load_synonyms(names_file):
    """Loads names.dmp and converts it into a pandas.DataFrame.
        Includes only names which are viewed as synonym by ncbi.
    """
    assert os.path.exists(names_file)
    # print("load synonyms")
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
    sci_df = df[df['name_class'] == 'synonym']
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


class Parser:
    """Reads in databases from ncbi to connect species names with the taxonomic identifier
    and the corresponding hierarchical information. It provides a much faster way to get those information then using
    web queries. We use those files to get independent from web requests to find those information (the implementation
    of it in BioPython was not really reliable).
    Nodes includes the hierarchical information, names the scientific names and ID's.
    The files need to be updated regularly, best way to always do it when a new blast database was loaded.
    """
    def __init__(self, names_file, nodes_file):
        self.names_file = names_file
        self.nodes_file = nodes_file
        self.initialize()

    def initialize(self):
        """ The data itself are not stored in __init__, as then the information will be pickled (which results in
        gigantic pickle file sizes).
        Instead every time the function is loaded after loading a pickle file, it will be 'initialized'.
        """
        global nodes
        nodes = load_nodes(self.nodes_file)
        global names
        names = load_names(self.names_file)
        global synonyms
        synonyms = load_synonyms(self.names_file)

    def get_rank(self, tax_id):
        """ Get rank for given ncbi tax id.
        """
        if nodes is None:
           self.initialize()
        rank = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        return rank


    def get_downtorank_id(self, tax_id, downtorank="species"):
        """ Recursive function to find the parent id of a taxon as defined by downtorank.
        """ 
        debug("get downtorank")
        if nodes is None:
            self.initialize()
        # following statement is to get id of taxa if taxa is higher ranked than specified
        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] != "species":
            if downtorank == "species":
                return tax_id
        if type(tax_id) != int:
            sys.stdout.write("WARNING: tax_id {} is no integer. Will convert value to int\n".format(tax_id))
            tax_id = int(tax_id)
        # debug(nodes[nodes["tax_id"] == tax_id])
        # debug(nodes[nodes['tax_id'] == tax_id]['rank'].values[0])
        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == downtorank:
            return tax_id
        elif nodes[nodes['tax_id'] == tax_id]['rank'].values[0] == "superkingdom":
            tax_id = 0
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

    def get_id_from_name(self, tax_name):
        """ Find the ID for a given taxonomic name.
        """
        if names is None:
            self.initialize()
        tax_name = tax_name.replace("_", " ")
        if len(tax_name.split(" ")) >= 2:
            if tax_name.split(" ")[1] == "sp.":
                tax_name = "{}".format(tax_name.split(" ")[0])
        try:
            tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
        except IndexError:
            if len(tax_name.split(" ")) == 3:
                tax_name = "{} {}-{}".format(tax_name.split(" ")[0], tax_name.split(" ")[1], tax_name.split(" ")[2])
                tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
            else:
                sys.stdout.write("Are you sure, its an accepted name and not a synonym?I look in the synonym table now")
                tax_id = self.get_id_from_synonym(tax_name)
        tax_id = int(tax_id)
        return tax_id

    def get_id_from_synonym(self, tax_name):
        """ Find the ID for a given taxonomic name, which is not an accepted name.
        """
        if names is None:
            self.initialize()
        tax_name = tax_name.replace("_", " ")
        try:
            tax_id = synonyms[synonyms["name_txt"] == tax_name]["tax_id"].values[0]
        except IndexError:
            if len(tax_name.split(" ")) == 3:
                tax_name = "{} {}-{}".format(tax_name.split(" ")[0], tax_name.split(" ")[1], tax_name.split(" ")[2])
                tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
            else:
                print("something else is going wrong: {}".format(tax_name))
        return tax_id

