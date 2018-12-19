import physcraper
import sys
from dendropy import DnaCharacterMatrix


def test_generate_ATT_from_phylesystem():
    seqaln = "tests/data/input.fas"
    study_id = "pg_873"
    tree_id = "tree1679"
    seqaln = "tests/data/minitest.fas"
    mattype = "fasta"
    workdir = "tests/output/opentree"
    configfi = "tests/data/remotencbi.config"


    sys.stdout.write("\nTesting 'generate_ATT_from_files (fromfile.py)'\n")

    conf = physcraper.ConfigObj(configfi, interactive=False)
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

    data_obj = physcraper.generate_ATT_from_phylesystem(aln=aln,
                                                    workdir=workdir,
                                                    config_obj=conf,
                                                    study_id=study_id,
                                                    tree_id=tree_id)

    data_obj == True
      