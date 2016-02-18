from physcraper import StudyInfo, PhyscraperSetup, PhyscraperScrape
import pickle
import sys
import os


study_id = "pg_873"
tree_id = "tree1679"
seqaln = "tests/data/minitest.fas"
mattype="fasta"
runname="fresh3"

info = StudyInfo(study_id, tree_id, seqaln, mattype)

sys.stdout.write("setting up StudyINfo\n")
sys.stdout.flush()
info = StudyInfo(study_id, tree_id, seqaln, mattype)
sys.stdout.write("setting up scrape instance\n")
sys.stdout.flush()

setup_pickfi = '{}/setup.p'.format(runname,runname)
scrape_pickfi = '{}/scrape.p'.format(runname,runname)

if os.path.isfile(scrape_pickfi):
    scrape = pickle.load(open(scrape_pickfi,'rb'))
elif os.path.isfile(setup_pickfi):
    setup = pickle.load(open(setup_pickfi,'rb'))
    scrape = PhyscraperScrape(setup)
else:
    setup = PhyscraperSetup(info, runname, configfi="config")
    setup.setup_physcraper()
    pickle.dump(setup, open(setup_pickfi, 'wb'))
    scrape = PhyscraperScrape(setup)

scrape.reset_markers()
sys.stdout.write("Instance set up\n")
sys.stdout.flush()
scrape.generate_streamed_alignment()
scrape.write_labelled()

