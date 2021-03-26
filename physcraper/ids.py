"""Link together NCBI and Open Tree identifiers and names, with Gen Bank information for updated sequences"""
import sys
import os
import time
from urllib.error import HTTPError
from Bio import Entrez
from physcraper import ncbi_data_parser, ConfigObj  # is the ncbi data parser class and associated functions
from physcraper.helpers import debug

_DEBUG = 0




class IdDicts():
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
            * self.ncbi_parser: for local blast,
               initializes the ncbi_parser class, that contains information about rank and identifiers
    """

    def __init__(self, configfile=None):
        """Generates a series of name disambiguation dicts"""
        if configfile is None:
            self.config = ConfigObj()
        elif isinstance(configfile, ConfigObj):
            self.config = configfile
        elif os.path.exists(configfile):
            self.config = ConfigObj(configfile)
        else:
            sys.stderr.write("Error reading config file {}\n".format(configfile))
            sys.exit()
        assert self.config
        self.ott_to_ncbi = {}
        self.ncbi_to_ott = {}  # used to get ott_id for new Genbank query taxa
        self.ott_to_name = {}  # used in add_otu to get name from otuId
        self.acc_ncbi_dict = {}  # filled by ncbi_parser (by subprocess in earlier versions of the code).
        self.spn_to_ncbiid = {}  # spn to ncbi_id, it's only fed by the ncbi_data_parser, but makes it faster
        self.ncbiid_to_spn = {}
        fi = open(self.config.ott_ncbi)
        # This is in the taxonomy folder of the repo, needs to be updated by devs when OpenTree taxonomy changes.
        # by downloading teh taxonomy and then using
        # grep ncbi: ottVERSION/taxonomy.tsv | sed -r -e
        # "s/([0-9]+).+?\|.+?\|(.+?)\|.+?\|.*ncbi:([0-9]+).*/\\1,\\3,\\2/" > physcraper/taxonomy/ott_ncbi
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
        else:  # ncbi parser contains information about spn, tax_id, and ranks
            debug("Config not remote {}".format(self.config.blast_loc))
            self.ncbi_parser = ncbi_data_parser.Parser(names_file=self.config.ncbi_names,
                                                       nodes_file=self.config.ncbi_nodes)
        self.acc_tax_seq_dict = {}
        self.full_seq_path = "{}/downloaded_ncbi_sequences".format(self.config.taxonomy_dir)


    def get_ncbiid_from_acc(self, acc):
        '''checks local dicts, and then runs eftech to get ncbi id for accession'''
        gb_id = acc
        if gb_id in self.acc_ncbi_dict:
            ncbi_id = self.acc_ncbi_dict[gb_id]
        elif gb_id in self.acc_tax_seq_dict:
            ncbi_id = self.acc_tax_seq_dict[gb_id]["^ncbi:taxon"]
        else:
            # Disabling unused variable bc not used locally but used by other functions later on
            taxid, taxname, seq = self.get_tax_seq_acc(gb_id) # pylint: disable=unused-variable
            ncbi_id = taxid
        return ncbi_id




 #removed function find_tax_id because it wasn't being used
    def get_tax_seq_acc(self, acc):
        """Pulls the taxon ID and the full sequences from NCBI"""
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
            ncbi_id = self.acc_tax_seq_dict[gb_id]["^ncbi:taxon"]
            seq = self.acc_tax_seq_dict[gb_id]["seq"]
        elif os.path.exists(seq_path):
            fi = open(seq_path)
            header = fi.readline().strip('>')
            assert header.split()[1].startswith('taxname:')
            tax_name = header.split()[1].strip('taxname:')
            ncbi_id = header.split()[2].strip('ncbi:')
            seq = "".join(fi.readlines())
        if seq is None:
            read_handle = self.entrez_efetch(gb_id)
            tax_name = ncbi_data_parser.get_ncbi_tax_name(read_handle)
            ncbi_id = ncbi_data_parser.get_ncbi_tax_id(read_handle)
            seq = read_handle[0][u'GBSeq_sequence']
            tax_name = tax_name.replace(" ", "_") #TODO check that searches are using names without spaces
            self.ncbiid_to_spn[ncbi_id] = tax_name
            self.acc_ncbi_dict[gb_id] = ncbi_id
            self.acc_tax_seq_dict[gb_id] = {'taxname':tax_name, "^ncbi:taxon":ncbi_id, 'seq':seq}#ToDo memory hog
            with open(seq_path, 'w') as fi:
                fi.write("> {} taxname:{} ncbi:{}\n".format(gb_id, tax_name, ncbi_id))
                fi.write(self.acc_tax_seq_dict[gb_id]['seq'])
        assert ncbi_id is not None
        return ncbi_id, tax_name, seq



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
                    raise e
            # break
        assert handle is not None, ("your handle file to access data from efetch does not exist. "
                                    "Likely an issue with the internet connection of ncbi. Try rerun...")
        read_handle = Entrez.read(handle)
        handle.close()
        return read_handle
