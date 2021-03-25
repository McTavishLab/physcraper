"""
Physcraper run Configuration object generator
"""
import sys
import os
import datetime
import configparser
import shutil
import wget


_DEBUG = 0




PHYSCRAPER_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
#sys.stdout.write(PHYSCRAPER_DIR)


def is_number(inputstr):
    """Test if string can be coerced to float"""
    try:
        float(inputstr)
        return True
    except ValueError:
        return False



class ConfigObj():
    """
    To build the class the following is needed:

      * **configfi**: a configuration file in a specific format, e.g. to read in self.e_value_thresh.

    During the initializing process the following self objects are generated:

      * **self.e_value_thresh**: the defined threshold for the e-value during Blast searches,
            check out: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
      * **self.hitlist_size**: the maximum number of sequences retrieved by a single blast search
      * **self.minlen**: value from 0 to 1. Defines how much shorter new seq can be compared to input
      * **self.trim_perc**: value that determines how many seq need to be present before the beginning
            and end of alignment will be trimmed
      * **self.maxlen**: max length for values to add to aln
      * **self.get_ncbi_taxonomy**: Path to sh file doing something...
      * **self.ott_ncbi**: file containing OTT id, ncbi and taxon name (??)
      * **self.email**: email address used for blast queries
      * **self.blast_loc**: defines which blasting method to use:

          * either web-query (=remote)
          * from a local blast database (=local)
      * **self.num_threads**: number of cores to be used during a run
      * **self.url_base**:

          * if blastloc == remote: it defines the url for the blast queries.
          * if blastloc == local: url_base = None
      * **self.delay**: defines when to reblast sequences in days
      * **optional self.objects**:

          * if blastloc == local:

              * self.blastdb: this defines the path to the local blast database
              * self.ncbi_nodes: path to 'nodes.dmp' file, that contains the hierarchical information
              * self.ncbi_names: path to 'names.dmp' file, that contains the different ID's
    """

    def __init__(self, configfile=None, run = True):
       # debug(configfi)
        if _DEBUG:
            sys.stdout.write("Building config object\n")
        self.run = run
        if configfile:
            self.set_defaults()
            self.read_config(configfile)
        else:
            sys.stdout.write("No config file, using defaults\n")
            self.set_defaults()
    def set_defaults(self):
        """ In the absence of an input configuration file, sets default values.
        """
        self.email = None
        self.e_value_thresh = 0.00001
        self.hitlist_size = 10
        self.blast_loc = 'remote'
        self.url_base = None
        self.blastdb = None
        self.num_threads = 4
        self.delay = 90
        self.spp_threshold = 5
        self.minlen = 0.8
        self.api_key = None
        self.maxlen = 1.2
        self.taxonomy_dir = ''
        self.check_taxonomy()
    def config_str(self):
        """ Write out the current config values.
        DOES NOT INCUDE SOME HIDDEN CONFIGUREABLE ATTRIBUTES
        """
        config_text = '''[blast]
Entrez.email = {email}
e_value_thresh = {e_val}
hitlist_size = {hls}
location = {bl}
localblastdb = {db}
num_threads = {nt}
delay = {delay}
[physcraper]
spp_threshold = {sppt}
min_length = {perc}
max_length = {maxlen}
'''.format(
    email=self.email,
    e_val=self.e_value_thresh,
    hls=self.hitlist_size,
    bl=self.blast_loc,
    db=self.blastdb,
    nt=self.num_threads,
    delay=self.delay,
    sppt=self.spp_threshold,
    perc=self.minlen,
    maxlen=self.maxlen)
        return config_text

    def check_taxonomy(self):
        """Locates a taxonomy directory in tthe phyysraper repo,
        or if not avail (often because module was pip installed),
        genertes one.
        """
        if os.path.exists(self.taxonomy_dir):
            pass
        elif os.path.exists("{}/taxonomy".format(PHYSCRAPER_DIR)):
            self.taxonomy_dir = "{}/taxonomy".format(PHYSCRAPER_DIR)
        else:
            if not os.path.exists('taxonomy'):
                os.mkdir("taxonomy")
            self.taxonomy_dir = os.path.abspath('taxonomy')
            sys.stdout.write("Using {} as taxonomy dir.".format(self.taxonomy_dir))
        self.ott_ncbi = "{}/ott_ncbi".format(self.taxonomy_dir)
        if not os.path.exists(self.ott_ncbi):
            sys.stdout.write(
                "downloading taxonomy from \
                https://raw.githubusercontent.com/McTavishLab/physcraper/main/taxonomy/ott_ncbi")
            wget.download(
                'https://raw.githubusercontent.com/McTavishLab/physcraper/main/taxonomy/ott_ncbi',
                out=self.taxonomy_dir)
        assert os.path.isfile(self.ott_ncbi), (
            "file `%s` does not exist" % self.ott_ncbi
            )

    def write_file(self, direc, filename="run.config"):
        """ writes config params to file
            * **direc**: path to write file
            * **filename**: filename to use. Default = run.config
        """
        config_text = self.config_str()
        fi = open("{}/{}".format(direc, filename), "w")
        fi.write(config_text)
        fi.close()

    def read_config(self, configfi):
        """ Reads configfile, and sets configuration params.
        any params not listed will be set to dafault values in set_default()
        * **configfile**: path to input file.
        """
        assert os.path.isfile(configfi), "file `%s` does not exist" % configfi
        config = configparser.ConfigParser()
        self.configfi = configfi
        config.read_file(open(configfi))

        # read in blast settings
        self.email = config["blast"].get("Entrez.email")
        #if not "@" in self.email:
            #sys.stderr.write(
            #    "your email `%s` does not have an @ sign. NCBI blast requests an email address.\n" % self.email)
        if config["blast"].get("Entrez.api_key"):
            self.api_key = config["blast"]["Entrez.api_key"]
            if self.api_key == 'None':
                self.api_key = None
        else:
            self.api_key = None
        self.e_value_thresh = config["blast"]["e_value_thresh"]
        assert is_number(self.e_value_thresh), (
            "value `%s` does not exists" % self.e_value_thresh
            )
        self.hitlist_size = int(config["blast"]["hitlist_size"])
        assert is_number(self.hitlist_size), (
            "value `%s`is not a number" % self.e_value_thresh
            )

        # read in settings for internal Physcraper processes
        if "taxonomy_path" in config["physcraper"].keys():
            self.taxonomy_dir = config["physcraper"]["taxonomy_path"]
        self.check_taxonomy()
        self.blast_loc = config["blast"]["location"]
        assert self.blast_loc in ["local", "remote"], (
            "your blast location `%s` is not remote or local" % self.email
        )
        if self.blast_loc == "local":
            self.blastdb = config["blast"]["localblastdb"]
            self.set_local()
        if self.blast_loc == "remote":
            self.url_base = config["blast"].get("url_base")
            if self.url_base == 'None':
                self.url_base = None
        if _DEBUG:
            sys.stdout.write("{}\n".format(self.email))
            #if self.blast_loc == "remote":
            #    sys.stdout.write("url base = {}\n".format(self.url_base))
            sys.stdout.write("{}\n".format(self.blast_loc))
            if self.blast_loc == "local":
                sys.stdout.write("local blast db {}\n".format(self.blastdb))
        self.num_threads = config["blast"].get("num_threads")
       # print("slurm threads")
       # print(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        if os.environ.get('SLURM_JOB_CPUS_PER_NODE'):
            self.num_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))
        self.delay = int(config["blast"]["delay"])
        assert is_number(self.delay), (
            "value `%s`is not a number" % self.delay
        )
        # #############
        # read in physcraper settings
        self.minlen = float(config["physcraper"]["min_length"])
        assert 0 < self.minlen <= 1, (
            "value `%s` is not between 0 and 1" % self.minlen
        )
        self.spp_threshold = int(config["physcraper"]["spp_threshold"])
        self.maxlen = float(config["physcraper"]["max_length"])
        assert self.maxlen > 1, (
            "value `%s` is not larger than 1" % self.maxlen
        )

    def set_local(self):
        """ Checks that all appropriate files etc are in place for local blast db.
        """
        if not self.run:
            return
        self.blast_loc = "local"
        self.ncbi_nodes = "{}/nodes.dmp".format(self.taxonomy_dir)
        self.ncbi_names = "{}/names.dmp".format(self.taxonomy_dir)
        assert(self.blastdb), "No blast db location set"
        if not os.path.isdir(self.blastdb):
            sys.stderr.write(
                "Local Blast DB not found at {},\
                please use a remote search, or update as described\
                in 'https://physcraper.readthedocs.io/en/latest/install.html#local-databases'\n".format(self.blastdb))
            sys.exit()
        if not os.path.exists("{}/nt.23.nhr".format(self.blastdb)):
            sys.stderr.write("Errors with local Blast DB at {}, \
                may be incomplete. please use a remote search, \
                or update as described in \
                'https://physcraper.readthedocs.io/en/latest/install.html#local-databases'\n".format(self.blastdb))
            sys.exit()
        else:
            download_date = os.path.getmtime("{}/nt.23.nhr".format(self.blastdb))
            download_date = datetime.datetime.fromtimestamp(download_date)
            today = datetime.datetime.now()
            time_passed = (today - download_date).days
            if time_passed >= 90:
                sys.stderr.write("Your databases might not be up to date anymore. \
                    You downloaded them {} days ago. \
                    Continuing, but perhaps use a remote search, \
                    or update as decribed in \
                    'https://physcraper.readthedocs.io/en/latest/install.html#local-databases'\n".format(time_passed))
        if not os.path.exists(self.ncbi_nodes):
            sys.stderr.write("NCBI taxonomy not found at {} - \
             To perform a local blast search, please update nodes and names.dmp, \
                as described in \
                'https://physcraper.readthedocs.io/en/latest/install.html#local-databases'\n".format(self.ncbi_nodes))

            sys.exit()
        else:
            download_date = os.path.getmtime(self.ncbi_nodes)
            download_date = datetime.datetime.fromtimestamp(download_date)
            today = datetime.datetime.now()
            time_passed = (today - download_date).days
            if time_passed >= 90:
                sys.stderr.write("Your taxonomy databases from NCBI were dowloaded {} days ago. \
                    To perform a local blast search, please update nodes and names.dmp, \
                    as described in \
                    'https://physcraper.readthedocs.io/en/latest/install.html#local-databases'\n".format(time_passed))
        assert(shutil.which("blastn")), "blastn  not found in path"
        self.url_base = None
