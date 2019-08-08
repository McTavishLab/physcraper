import sys
import re
import os
import subprocess
import datetime
import glob
import json
import configparser
import pickle
import random
import time
import csv

from physcraper import ncbi_data_parser  # is the ncbi data parser class and associated functions

from physcraper.helpers import cd, debug

_DEBUG = 0


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


def get_user_input():
    """Asks for yes or no user input.

    :return: user input
    """
    debug("get user input")
    is_valid = 0
    x = None
    while not is_valid:
        try:
            x = raw_input("Please write either 'yes' or 'no': ")
            if x == 'yes' or x == 'no':
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
            else:
                print("'%s' is not a valid answer." % x)
        except ValueError as e:
            print("'%s' is not a valid answer." % e.args[0].split(": ")[1])
    return x


class ConfigObj(object):
    """
    To build the class the following is needed:

      * **configfi**: a configuration file in a specific format, e.g. to read in self.e_value_thresh.
                
        The file needs to have a heading of the format: [blast] and then somewhere below that heading a string e_value_thresh = value

      * **interactive**: defaults to True, is used to interactively update the local blast databases

    During the initializing process the following self objects are generated:

      * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches, check out: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
      * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
      * **self.seq_len_perc**: value from 0 to 1. Defines how much shorter new seq can be compared to input
      * **self.trim_perc**: value that determines how many seq need to be present before the beginning and end of alignment will be trimmed
      * **self.maxlen**: max length for values to add to aln
      * **self.get_ncbi_taxonomy**: Path to sh file doing something...
      * **self.phylesystem_loc**: defines which phylesystem for OpenTree datastore is used. The default is api, but can run on local version too.
      * **self.ott_ncbi**: file containing OTT id, ncbi and taxon name (??)
      * **self.id_pickle**: path to pickle file
      * **self.email**: email address used for blast queries
      * **self.blast_loc**: defines which blasting method to use:

          * either web-query (=remote)
          * from a local blast database (=local)
      * **self.num_threads**: number of cores to be used during a run
      * **self.url_base**: 

          * if blastloc == remote: it defines the url for the blast queries.
          * if blastloc == local: url_base = None
      * **self.unmapped**: used for OToL original tips that can not be assigned to a taxon

          * keep: keep the unmapped taxa and asign them to life
          * remove: remove the unmapped taxa from aln and tre
      * **self.delay**: defines when to reblast sequences in days
      * **optional self.objects**:

          * if blastloc == local:

              * self.blastdb: this defines the path to the local blast database
              * self.ncbi_parser_nodes_fn: path to 'nodes.dmp' file, that contains the hierarchical information
              * self.ncbi_parser_names_fn: path to 'names.dmp' file, that contains the different ID's
    """

    def __init__(self, configfi, interactive=None):
       # debug(configfi)
        assert os.path.isfile(configfi)
        if _DEBUG:
            sys.stdout.write("Building config object\n")
        assert os.path.isfile(configfi), "file `%s` does not exists" % configfi
        config = configparser.ConfigParser()
        config.read_file(open(configfi))
        
        # read in blast settings
        self.email = config["blast"]["Entrez.email"]
        assert "@" in self.email, "your email `%s` does not have an @ sign" % self.email
        
        self.e_value_thresh = config["blast"]["e_value_thresh"]
        assert is_number(self.e_value_thresh), (
                "value `%s` does not exists" % self.e_value_thresh
        )
        self.hitlist_size = int(config["blast"]["hitlist_size"])
        assert is_number(self.hitlist_size), (
                "value `%s`is not a number" % self.e_value_thresh
        )
        self.blast_loc = config["blast"]["location"]
        assert self.blast_loc in ["local", "remote"], (
                "your blast location `%s` is not remote or local" % self.email
        )        
        if self.blast_loc == "local":
            self.blastdb = config["blast"]["localblastdb"]
            self.url_base = None
            self.ncbi_parser_nodes_fn = config["ncbi_parser"]["nodes_fn"]
            self.ncbi_parser_names_fn = config["ncbi_parser"]["names_fn"]
        if self.blast_loc == "remote":
            self.url_base = config["blast"].get("url_base")
        if _DEBUG:
            sys.stdout.write("{}\n".format(self.email))
            if self.blast_loc == "remote":
                sys.stdout.write("url base = {}\n".format(self.url_base))
            sys.stdout.write("{}\n".format(self.blast_loc))
            if self.blast_loc == "local":
                sys.stdout.write("local blast db {}\n".format(self.blastdb))
        self.num_threads = config["blast"].get("num_threads")
       # print("slurm threads")
       # print(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        if os.environ.get('SLURM_JOB_CPUS_PER_NODE'):
            self.num_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))

      #  debug(self.num_threads)
        self.gb_id_filename = config["blast"].get("gb_id_filename", False)
        if self.gb_id_filename is not False:
            if self.gb_id_filename == "True" or self.gb_id_filename == "true":
                self.gb_id_filename = True
            else:
                self.gb_id_filename = False
     #   debug("shared blast folder? {}".format(self.gb_id_filename))
        self.delay = int(config["blast"]["delay"])
        assert is_number(self.delay), (
                "value `%s`is not a number" % self.delay
        )       
        # #############
        # read in physcraper settings       
        self.unmapped = config["physcraper"]["unmapped"]
        assert self.unmapped in ["remove", "keep"], (
                "your unmapped statement `%s` in the config file is not remove or keep"
                % self.unmapped
        )
        self.seq_len_perc = float(config["physcraper"]["seq_len_perc"])
        assert 0 < self.seq_len_perc <= 1, (
                "value `%s` is not between 0 and 1" % self.seq_len_perc
        )
        self.trim_perc = float(config["physcraper"]["trim_perc"])
        assert 0 < self.trim_perc < 1, (
                "value `%s` is not between 0 and 1" % self.trim_perc
        )
        self.maxlen = float(config["physcraper"]["max_len"])
        assert 1 < self.maxlen, (
                "value `%s` is not larger than 1" % self.maxlen
        )
        

        # read in settings for internal Physcraper processes
        self.phylesystem_loc = config["phylesystem"]["location"]
        assert self.phylesystem_loc in ["local", "api"], \
            (
                "phylesystem location must be either local or api")  # default is api, but can run on local version of OpenTree datastore
        self.ott_ncbi = config["taxonomy"][
            "ott_ncbi"
        ]
        assert os.path.isfile(self.ott_ncbi), (
                "file `%s` does not exists" % self.ott_ncbi
        )
        # rewrites relative path to absolute path so that it behaves when changing dirs
        self.id_pickle = os.path.abspath(config["taxonomy"]["id_pickle"])
        
        ####
        # check database status
        if interactive is None:
            interactive = sys.stdin.isatty()
            if interactive is False:
                sys.stdout.write("REMEMBER TO UPDATE THE NCBI DATABASES REGULARLY!!\n")
        if interactive is True:
            self._download_ncbi_parser()
            self._download_localblastdb()
#        debug("check db file status?: {}".format(interactive))

    def _download_localblastdb(self):
        """Check if files are present and if they are uptodate.
        If not files will be downloaded.
        """
        if self.blast_loc == "local":
            # next line of codes exists to have interactive mode enabled while testing
            # this allows to not actually have a local ncbi database downloaded
            if not os.path.isfile("{}/empty_local_db_for_testing.nhr".format(self.blastdb)):
                if not os.path.isfile("{}/nt.69.nhr".format(self.blastdb)):
                    print("To run local blast queries, download the blast data basefrom ncbi. See https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/")
                    self.blast_loc = "remote"
                else:
                    download_date = os.path.getmtime("{}/nt.60.nhr".format(self.blastdb))
                    download_date = datetime.datetime.fromtimestamp(download_date)
                    today = datetime.datetime.now()
                    time_passed = (today - download_date).days
                    if time_passed >= 90:
                        print("Your databases might not be uptodate anymore. You downloaded them {} days ago. "
                              "Do you want to update the blast databases from ncbi? Note: This is a US government website! "
                              "You agree to their terms".format(time_passed))
                        x = get_user_input()
                        if x == "yes":
                            with cd(self.blastdb):
                                os.system("update_blastdb nt")
                                os.system("cat *.tar.gz | tar -xvzf - -i")
                                os.system("update_blastdb taxdb")
                                os.system("gunzip -cd taxdb.tar.gz | (tar xvf - )")
                                os.system("rm *.tar.gz*")
                        elif x == "no":
                            print("You did not agree to update data from ncbi. Old database files will be used.")
                        else:
                            print("You did not type 'yes' or 'no'!")

    def _download_ncbi_parser(self):
        """Check if files are present and if they are up to date.
        If not files will be downloaded. 
        """
        if self.blast_loc == "local":
            if not os.path.isfile(self.ncbi_parser_nodes_fn):
                print("Do you want to download taxonomy databases from ncbi? Note: This is a US government website! "
                      "You agree to their terms")
                x = get_user_input()
                if x == "yes":
                    os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./tests/data/")
                    os.system("gunzip -f -cd ./taxonomy/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                    os.system("mv nodes.dmp ./taxonomy")
                    os.system("mv names.dmp ./taxonomy/")
                    os.system("rm taxdump.tar.gz")
                elif x == "no":
                    print("You did not agree to download data from ncbi. Program will default to blast web-queries.")
                    print("This is slow and crashes regularly!")
                    self.blast_loc = "remote"
                else:
                    print("You did not type yes or no!")
            else:
                download_date = os.path.getmtime(self.ncbi_parser_nodes_fn)
                download_date = datetime.datetime.fromtimestamp(download_date)
                today = datetime.datetime.now()
                time_passed = (today - download_date).days
                if time_passed >= 90:
                    print("Do you want to update taxonomy databases from ncbi? Note: This is a US government website! "
                          "You agree to their terms")
                    x = get_user_input()
                    if x == "yes":
                        os.system("wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' -P ./taxonomy/")
                        os.system("gunzip -f -cd ./tests/data/taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)")
                        os.system("mv nodes.dmp ./taxonomy/")
                        os.system("mv names.dmp ./taxonomy/")
                    elif x == "no":
                        print("You did not agree to update data from ncbi. Old database files will be used.")
                    else:
                        print("You did not type yes or no!")
