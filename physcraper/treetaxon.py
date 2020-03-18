import sys
import re
import json
import os
from dendropy import Tree
#from opentree_helpers import bulk_tnrs_load





class TreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match.

    """
    def __init__(self, otu_json, treefrom = 'synth', 
                 schema='newick', taxon_namespace=None):
        if treefrom == 'synth':
            pass
        elif os.path.exists(treefrom):
            self.tre = Tree.get(path=treefrom,
                                schema=schema,
                                preserve_underscores=True)
        if isinstance(otu_json, dict):
            self.otu_dict = otu_dict
        elif isinstance(otu_json, str):
            assert os.path.exists(otu_json)
            with open(otu_json) as data_file:
                input_dict = json.load(data_file)
                if input_dict.keys() == [u'mappingHints', u'names', u'metadata']:
                    self.otu_dict = bulk_tnrs_load(otu_json)
                else:
                    self.otu_dict = input_dict
        self._reconcile_names()
    def _reconcile_names(self):
        """It rewrites some tip names, which kept being an issue when it starts with a number at the beginning.
        Then somehow a n was added to the tip names.

        :return: replaced tip names
        """
        i = 1
        for tax in self.tre.taxon_namespace:
            if tax.label in self.otu_dict.keys():
                pass
            else:
                found_label = 0
                match = re.match("'n[0-9]{1,3}", tax.label)
                newname = ""
                if match:
                    newname = tax.label[2:]
                    newname = newname[:-1]
                for otu in self.otu_dict:
                    original = self.otu_dict[otu].get("^ot:originalLabel")
                    if original == tax.label or original == newname:
                        tax.label = otu
                        found_label = 1
                if found_label == 0:
                    sys.stderr.write("could not match tiplabel {} or {} to an OTU\n".format(tax.label, newname))
                    otulab = "otu"+tax.label
                    self.otu_dict[otulab]["^ot:originalLabel"] = tax.label
                    tax.label = otu_lab
        for tax in self.tre.taxon_namespace:
              assert tax.label in self.otu_dict          
    def write_labelled(self, filepath,  label='^ot:ottTaxonName', norepeats=True, add_gb_id=False):
        """output tree and alignment with human readable labels
        Jumps through a bunch of hoops to make labels unique.

        :param label: which information shall be displayed in labelled files: possible options:
                    '^ot:ottTaxonName', '^user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"
        :param treepath: optional: full file name (including path) for phylogeny
        :param norepeats: optional: if there shall be no duplicate names in the labelled output files
        :param add_gb_id: optional, to supplement tiplabel with corresponding GenBank sequence identifier
        :return: writes out labelled phylogeny and alignment to file
        """
        #debug("write labelled files")
        assert label in ['^ot:ottTaxonName', '^user:TaxonName', '^physcraper:TaxonName',
                         "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon", '^ncbi:accession']
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:
            new_label = self.otu_dict[taxon.label].get(label, None)
            if new_label is None:
                if self.otu_dict[taxon.label].get("^ot:originalLabel"):
                    new_label = "orig_{}".format(self.otu_dict[taxon.label]["^ot:originalLabel"])
                else:
                    new_label = "ncbi_{}_ottname_{}".format(self.otu_dict[taxon.label].get("^ncbi:taxon", "unk"),
                                                            self.otu_dict[taxon.label].get('^physcraper:TaxonName', "unk"))
                new_label = str(new_label).replace(' ', '_')
            else:
                if new_label in new_names and norepeats:
                    new_label = "_".join([new_label, taxon.label])
            taxon.label = new_label
            new_names.add(new_label)
        tmp_tre.write(path=filepath,
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)
