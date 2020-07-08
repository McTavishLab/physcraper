import sys
import os
import json
import time


if sys.version_info < (3,):
    from urllib2 import HTTPError
else:
    from urllib.error import HTTPError

from Bio import Entrez

from physcraper import ncbi_data_parser, ConfigObj  # is the ncbi data parser class and associated functions
from physcraper.helpers import debug

_DEBUG = 0




class IdDicts(object):
    """Class contains different taxonomic identifiers and helps to find the corresponding ids between ncbi and OToL

        To build the class the following is needed:

          * **config_obj**: Object of class config (see above)
          * **workdir**: the path to the assigned working directory

        During the initializing process the following self objects are generated:

          * **self.workdir**: contains path of working directory
          * **self.config**: contains the Config class object
          * **self.ott_to_ncbi**: dictionary

              * key: OToL taxon identifier
              * value: ncbi taxon identifier
          * **self.ncbi_to_ott**: dictionary

              * key: OToL taxon identifier
              * value: ncbi taxon identifier
          * **self.ott_to_name**: dictionary

              * key: OToL taxon identifier
              * value: OToL taxon name
          * **self.acc_ncbi_dict**: dictionary

              * key: Genbank identifier
              * value: ncbi taxon identifier
          * **self.spn_to_ncbiid**: dictionary

              * key: OToL taxon name
              * value: ncbi taxon identifier
          * **self.ncbiid_to_spn**: dictionary

              * key: ncbi taxon identifier
              * value: ncbi taxon name

          user defined list of mrca OTT-ID's #TODO this is flipped form the dat aobj .ott_mrca. On purpose?
         #reomved mrca's from ida, and put them into scrape object

          * **Optional**:

              * depending on blasting method:
               * self.ncbi_parser: for local blast, initializes the ncbi_parser class, that contains information about rank and identifiers
               * self.otu_rank: for remote blast to store the rank information
    """

    def __init__(self, configfile = None, workdir=None):
        """Generates a series of name disambiguation dicts"""
        if configfile == None:
            self.config = ConfigObj()
        elif isinstance(configfile, ConfigObj):
            self.config = configfile
        elif os.path.exists(configfile):
            self.config = ConfigObj(configfile)
        else:
            sys.stderr.write("Error reading config file\n".format(configfile))
            sys.exit()
        assert self.config
        self.ott_to_ncbi = {} 
        self.ncbi_to_ott = {}  # used to get ott_id for new Genbank query taxa
        self.ott_to_name = {}  # used in add_otu to get name from otuId
        self.acc_ncbi_dict = {}  # filled by ncbi_parser (by subprocess in earlier versions of the code).
        self.spn_to_ncbiid = {}  # spn to ncbi_id, it's only fed by the ncbi_data_parser, but makes it faster
        self.ncbiid_to_spn = {} #TODO when is this generated? MK: well, here. it is filled with information from genbank to speed up translation between ncbi_taxon_ids and names. similar to  acc_ncbi_dict and spn_to_ncbiid.
        tax_folder = self.config.taxonomy_dir
        fi = open(self.config.ott_ncbi)  # This is in the taxonomy folder of the repo, needs to be updated by devs when OpenTree taxonomy changes.
        for lin in fi:  
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
            self.ott_to_name[int(lii[0])] = lii[2].strip()  # todo merge into ott_to_ncbi?
        fi.close()
        assert len(self.ott_to_ncbi) > 0
        assert len(self.ott_to_name) > 0
        assert len(self.ncbi_to_ott) > 1000
        if self.config.blast_loc == 'remote':
            debug("Config remote {}".format(self.config.blast_loc))
            self.otu_rank = {}  # used only for web queries - contains taxonomic hierarchy information
        else:  # ncbi parser contains information about spn, tax_id, and ranks
            debug("Config not remote {}".format(self.config.blast_loc))
            self.ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_names,
                                                       nodes_file=self.config.ncbi_nodes)
        self.acc_tax_seq_dict = {}
        self.full_seq_path = "{}/downloaded_ncbi_sequences".format(self.config.taxonomy_dir)


    def get_ncbiid_from_acc(self, acc):
        '''checks local dicts, and then runs eftech to get ncbi id for accession'''
        gb_id = acc
        if gb_id in self.acc_ncbi_dict:#TODO if the accession number and tax id are here, does that mean the name is in ncbiid_to_spn?
            ncbi_id = self.acc_ncbi_dict[gb_id]
        elif gb_id in self.acc_tax_seq_dict:
            ncbi_id = self.acc_tax_seq_dict[gb_id]["^ncbi:taxon"]
        else:
            taxid, taxname, seq = self.get_tax_seq_acc(gb_id)
            ncbi_id = taxid
        return ncbi_id




 #removed function find_tax_id because it wasn't being used
    def get_tax_seq_acc(self, acc):
        if not os.path.exists(self.full_seq_path):
            os.makedirs(self.full_seq_path)
        gb_id = acc
        if len(gb_id.split(".")) == 1:
            debug("accession number {} not recognized".format(gb_id))
            return None, None, None
        seq_path = "{}/{}_remote.fasta".format(self.full_seq_path, gb_id)
        ncbi_id = tax_name = seq = None
        if gb_id in self.acc_tax_seq_dict:
            tax_name = self.acc_tax_seq_dict[gb_id]["taxname"]
            ncbi_id = self.acc_tax_seq_dict[gb_id]["^ncbi:"]
            seq = self.acc_tax_seq_dict[gb_id]["seq"]
        elif os.path.exists(seq_path):
            fi = open(seq_path)
            header = fi.readline().strip('>')
            #try:
            assert(header.split()[1].startswith('taxname:'))
            tax_name = header.split()[1].strip('taxname:')
            ncbi_id = header.split()[2].strip('ncbi:')
            seq = "".join(fi.readlines())
#            except IndexError:
 #               print("IndexError")
  #              pass
        if seq == None:
            read_handle = self.entrez_efetch(gb_id)
            tax_name = ncbi_data_parser.get_ncbi_tax_name(read_handle)
            ncbi_id =  ncbi_data_parser.get_ncbi_tax_id(read_handle)
            seq = read_handle[0][u'GBSeq_sequence']
            tax_name = tax_name.replace(" ","_") #TODO check that searches are using names without spaces 
            self.ncbiid_to_spn[ncbi_id] = tax_name
            self.acc_ncbi_dict[gb_id] = ncbi_id
            self.acc_tax_seq_dict[gb_id] = {'taxname':tax_name, "^ncbi:taxon":ncbi_id, 'seq':seq} #This is going to be a memory hog...
            with open(seq_path, 'w') as fi:
                fi.write("> {} taxname:{} ncbi:{}\n".format(gb_id, tax_name, ncbi_id))
                fi.write(self.acc_tax_seq_dict[gb_id]['seq'])
        assert ncbi_id is not None
        return ncbi_id, tax_name, seq




    def find_name_otu(self, otu_dict_entry=None):
        """ Find the taxon name in the  otu_dict entry or of a Genbank accession number.
        If not already known it will ask ncbi using the accession number

        :param otu_dict_entry: otu_dict entry
        :param acc: Genbank accession number
        :return: ncbi taxon name
        """
        # debug("find_name")
        inputinfo = False
        if otu_dict_entry is not None:
            inputinfo = True
        assert inputinfo is True
        tax_name = None
        ncbi_id = None
        if otu_dict_entry:
            # debug(otu_dict_entry)
            if "^physcraper:TaxonName" in otu_dict_entry:
                tax_name = otu_dict_entry["^physcraper:TaxonName"]
            elif "^ot:ottTaxonName" in otu_dict_entry:
                tax_name = otu_dict_entry["^ot:ottTaxonName"]
            elif "^user:TaxonName" in otu_dict_entry:
                tax_name = otu_dict_entry["^user:TaxonName"]
        assert tax_name is not None
        tax_name = tax_name.replace(" ", "_")
        return tax_name


    def entrez_efetch(self, gb_id):
        """ Wrapper function around efetch from ncbi to get taxonomic information if everything else is failing.
            Also used when the local blast files have redundant information to access the taxon info of those sequences.
        It adds information to various id_dicts.

        :param gb_id: Genbank identifier
        :return: read_handle
        """
        tries = 10
        Entrez.email = self.config.email
        if self.config.api_key:
            Entrez.api_key = self.config.api_key
        handle = None

        # method needs delay because of ncbi settings
        for i in range(tries):
            try:
                # print("try")
                delay = 1.0
                previous = time.time()
                while True:
                    current = time.time()
                    wait = previous + delay - current
                    if wait > 0:
                        # print("if", wait)
                        time.sleep(wait)
                        previous = current + wait
                    else:
                        # print("else", wait)
                        previous = current
                    if delay + .5 * delay <= 5:
                        # print("if2", delay)
                        delay += .5 * delay
                    else:
                        # print("else2",  delay)
                        delay = 5
                    # print("read handle")
                    handle = Entrez.efetch(db="nucleotide", id=gb_id, retmode="xml")
                    assert handle is not None, ("your handle file to access data from efetch does not exist. "
                                                "Likely an issue with the internet connection of ncbi. Try rerun...")
                    read_handle = Entrez.read(handle)
                    handle.close()

                    return read_handle
            except (IndexError, HTTPError) as e:
                if i < tries - 1:  # i is zero indexed
                    continue
                else:
                    raise
            # break
        assert handle is not None, ("your handle file to access data from efetch does not exist. "
                                    "Likely an issue with the internet connection of ncbi. Try rerun...")
        read_handle = Entrez.read(handle)
        handle.close()
        return read_handle


    def dump(self, filename=None):
        if filename:
            ofi = open(filename, "wb")
        else:
            ofi = open("id_pickle.p", "wb")
        pickle.dump(self, ofi)

