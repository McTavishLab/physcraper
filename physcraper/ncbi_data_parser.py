"""uses ncbi databases to easily retrieve taxonomic information.

parts are altered from https://github.com/zyxue/ncbitax2lin/blob/master/ncbitax2lin.py
"""

import os
import sys
import pandas as pd
from physcraper.helpers import debug



_DEBUG = 0


def get_acc_from_blast(query_string):
    """
    Get the accession number from a blast query.
    :param query_string: string that contains acc and gi from local blast query result
    :return: gb_acc

    """
    if len(query_string.split("|")) >= 3:
        gb_acc = query_string.split("|")[3]
    else:
        gb_acc = query_string.split("|")[0]
    if len(gb_acc.split(".")) < 2:
        sys.stderr.write("query string {} does not contain a Genbank accession number. Skipping\n".format(query_string))
        return None
    assert len(gb_acc.split(".")) >= 2, (len(gb_acc.split(".")), gb_acc)
    return gb_acc

def get_gi_from_blast(query_string):
    """
    Get the gi number from a blast query.
    Get acc is more difficult now, as new seqs not always have gi number, then query changes.

    If not available return None.

    :param query_string: string that contains acc and gi from local blast query result
    :return: gb_id if available
    """
    if len(query_string.split("|")) >= 3:
        gb_id = query_string.split("|")[1]
        assert len(gb_gi.split(".")) < 2, (len(gb_gi.split(".")), gb_gi)
        assert gb_gi.isdigit() is True
        return gb_gi
     else:
        return None

def get_tax_info_from_acc(gb_id, ids_obj):
    '''takes an accession number and returns the ncbi_id and the taxon name'''
#    debug("Getting tax info from acc {}".format(gb_id))
    ncbi_id = ids_obj.get_ncbiid_from_acc(gb_id)
    tax_name = ids_obj.ncbiid_to_spn.get(ncbi_id)
    if ncbi_id is None:
        sys.stderr.write("Failed to get information for sequence with accession number {}\n".format(gb_id))
    return ncbi_id, tax_name


def get_ncbi_tax_id(handle):
    """Get the taxon ID from ncbi. ONly used for web queries

    :param handle: NCBI read.handle
    :return: ncbi_id
    """
    ncbi_id = None
    gb_list = handle[0]["GBSeq_feature-table"][0]["GBFeature_quals"]
    for item in gb_list:
        if item[u"GBQualifier_name"] == "db_xref":
            if item[u"GBQualifier_value"][:5] == "taxon":
                ncbi_id = int(item[u"GBQualifier_value"][6:])
                break
            else:
                continue
    return ncbi_id


def get_ncbi_tax_name(handle):
    """Get the sp name from ncbi.
    Could be replaced by direct lookup to ott_ncbi.

    :param handle: NCBI read.handle
    :return: ncbi_spn
    """
    ncbi_sp = None
    gb_list = handle[0]["GBSeq_feature-table"][0]["GBFeature_quals"]
    for item in gb_list:
        if item[u"GBQualifier_name"] == "organism":
            ncbi_sp = str(item[u"GBQualifier_value"])
            ncbi_sp = ncbi_sp.replace(" ", "_")
    return ncbi_sp


def strip(input):
    """ Strips of blank characters from string in pd dataframe.
    """
    if isinstance(input, str):
        return input.strip()
    else:
        return input


def load_nodes(nodes_file):
    """ Loads nodes.dmp and converts it into a pandas.DataFrame.
    Contains the information about the taxonomic hierarchy of names.
    """
    # print(nodes_file)
    assert os.path.exists(nodes_file), (
        "file `%s` does not exist. Make sure you downloaded the "
        "databases from ncbi." % nodes_file
    )
    df = pd.read_csv(
        nodes_file,
        sep="|",
        header=None,
        index_col=False,
        names=[
            "tax_id",
            "parent_tax_id",
            "rank",
            "embl_code",
            "division_id",
            "inherited_div_flag",
            "genetic_code_id",
            "inherited_GC__flag",
            "mitochondrial_genetic_code_id",
            "inherited_MGC_flag",
            "GenBank_hidden_flag",
            "hidden_subtree_root_flag",
            "comments",
        ],
    )
    # To get rid of flanking tab characters
    df["rank"] = df["rank"].apply(strip)
    df["embl_code"] = df["embl_code"].apply(strip)
    df["comments"] = df["comments"].apply(strip)
    return df


def load_names(names_file):
    """ Loads names.dmp and converts it into a pandas.DataFrame.
    Includes only names which are accepted as scientific name by ncbi.
    """
    assert os.path.exists(names_file), (
        "file `%s` does not exist. Make sure you downloaded the "
        "databases from ncbi." % names_file
    )
    df = pd.read_csv(
        names_file,
        sep="|",
        header=None,
        index_col=False,
        names=["tax_id", "name_txt", "unique_name", "name_class"],
    )
    df["name_txt"] = df["name_txt"].apply(strip)
    df["unique_name"] = df["unique_name"].apply(strip)
    df["name_class"] = df["name_class"].apply(strip)
    sci_df = df[df["name_class"] == "scientific name"]
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


def load_synonyms(names_file):
    """Loads names.dmp and converts it into a pandas.DataFrame.
        Includes only names which are viewed as synonym by ncbi.
    """
    assert os.path.exists(names_file), (
        "file `%s` does not exist. Make sure you downloaded "
        "the databases from ncbi." % names_file
    )
    # print("load synonyms")
    df = pd.read_csv(
        names_file,
        sep="|",
        header=None,
        index_col=False,
        names=["tax_id", "name_txt", "unique_name", "name_class"],
    )
    df["name_txt"] = df["name_txt"].apply(strip)
    df["unique_name"] = df["unique_name"].apply(strip)
    df["name_class"] = df["name_class"].apply(strip)
    sci_df = df[df["name_class"] == "synonym"]
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
        # self.initialize()
        self.nodes = load_nodes(self.nodes_file)
        self.names = load_names(self.names_file)
        self.synonyms = load_synonyms(self.names_file)

    def get_rank(self, tax_id):
        """ Get rank for given ncbi tax id.
        """
        if tax_id == None:
            rank = "unassigned"
        else:
            nodes = self.nodes
            rank = nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
        return rank

    def get_downtorank_id(self, tax_id, downtorank="species"):
        """ Recursive function to find the parent id of a taxon as defined by downtorank.
        """
        nodes = self.nodes
        if not isinstance(tax_id, int):
#        debug("get downtorank")
           # sys.stdout.write(
           #     "WARNING: tax_id {} is no integer. Will convert value to int\n".format(
           #         tax_id
           #     )
           # )
            tax_id = int(tax_id)
        # debug(downtorank)
        # following statement is to get id of taxa if taxa is higher ranked than specified
        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] != "species":
            if downtorank == "species":
                if (nodes[nodes["tax_id"] == tax_id]["rank"].values[0] != "varietas"
                        and nodes[nodes["tax_id"] == tax_id]["rank"].values[0]
                        != "subspecies"):
                    return tax_id
        if nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == downtorank:
            # debug("found right rank")
            return tax_id
        elif nodes[nodes["tax_id"] == tax_id]["rank"].values[0] == "superkingdom":
            tax_id = 0
            return tax_id
        else:
            parent_id = int(nodes[nodes["tax_id"] == tax_id]["parent_tax_id"].values[0])
            return self.get_downtorank_id(parent_id, downtorank)

    def match_id_to_mrca(self, tax_id, mrca_id):
        """ Recursive function to find out if tax_id is part of mrca_id.
        """
        # debug("match_id_to_mrca")
        nodes = self.nodes
       # debug("testing if {} within {}".format(tax_id, mrca_id))
        current_id = int(tax_id)
        mrca_id = int(mrca_id)
        #debug([rank_mrca_id, rank_tax_id])
        while current_id:
            if current_id == mrca_id:
                # debug("found right rank")
                return True
            elif current_id == 1:
   #             debug("current id is: {}".format(current_id))
                return False
            elif current_id == 0:
                debug("current id is: {}, in search for {} in {}".format(current_id, tax_id, mrca_id))
                return False
            else: #try parent
                try:
                    current_id = int(nodes[nodes["tax_id"] == current_id]["parent_tax_id"].values[0])
                    # to get the except:
                    # current_id = 17043521
                    # nodes = ncbitax.nodes  ## ncbitax object from test_ncbi_parser.py
                except IndexError:
                    sys.stderr.write("no parent found for ncbi:id {}".format(current_id))
                    return False
#                debug("parent id is: {}".format(current_id))

    def get_name_from_id(self, tax_id):
        """ Find the scientific name for a given ID.
        """
        try:
            names = self.names
            if tax_id == 0:
                tax_name = "unidentified"
            else:
                tax_name = names[names["tax_id"] == tax_id]["name_txt"]
                tax_name = tax_name.values[0].replace(" ", "_")
                tax_name = tax_name.strip()
        except IndexError:
            sys.stdout.write("tax_id {} unknown by ncbi_parser files (names.dmp)\n".format(tax_id))
            tax_name = "unknown_{}".format(tax_id)
            if os.path.exists("ncbi_id_unknown.err"):
                fn = open("ncbi_id_unknown.err", "a")
                fn.write("{}".format(tax_id))
                fn.close()
            else:
                fn = open("ncbi_id_unknown.err", "w")
                fn.write("{}".format(tax_id))
                fn.close()
        return tax_name

    def get_id_from_name(self, tax_name):
        """ Find the ID for a given taxonomic name.
        """
        names = self.names
        org_tax = tax_name
        tax_name = tax_name.replace("_", " ")
        if len(tax_name.split(" ")) >= 2:
            if tax_name.split(" ")[1] == "sp.":
                tax_name = "{}".format(tax_name.split(" ")[0])
        try:
            tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
        except IndexError:
            if len(tax_name.split(" ")) == 3:
                tax_name = "{} {}-{}".format(
                    tax_name.split(" ")[0],
                    tax_name.split(" ")[1],
                    tax_name.split(" ")[2],
                )
                tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
                sys.stdout.write(
                    "tax_name {} unknown, modified to {} worked.\n".format(org_tax, tax_name)
                )
            else:
                sys.stdout.write(
                    "Are you sure, its an accepted name and not a synonym: {}? "
                    "I look in the synonym table now.\n".format(tax_name)
                )
                tax_id = self.get_id_from_synonym(tax_name)
        tax_id = int(tax_id)
        return tax_id

    def get_id_from_synonym(self, tax_name):
        """ Find the ID for a given taxonomic name, which is not an accepted name.
        """
        names = self.names
        synonyms = self.synonyms
        tax_name = tax_name.replace("_", " ")
        try:
            tax_id = synonyms[synonyms["name_txt"] == tax_name]["tax_id"].values[0]
        except IndexError:
            if len(tax_name.split(" ")) == 3:
                tax_name = "{} {}-{}".format(
                    tax_name.split(" ")[0],
                    tax_name.split(" ")[1],
                    tax_name.split(" ")[2],
                )
                tax_id = names[names["name_txt"] == tax_name]["tax_id"].values[0]
            else:
                sys.stderr.write("ncbi taxon name unknown by parser files: {}, taxid set to 0.\n".format(tax_name))
                tax_id = 0
                if os.path.exists("ncbi_name_unknown.err"):
                    fn = open("ncbi_name_unknown.err", "a")
                    fn.write("{}".format(tax_id))
                    fn.close()
                else:
                    fn = open("ncbi_name_unknown.err", "w")
                    fn.write("{}".format(tax_id))
                    fn.close()
        return tax_id
