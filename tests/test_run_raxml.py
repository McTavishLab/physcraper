import pickle
import sys
import os
from physcraper import ConfigObj, PhyscraperScrape, IdDicts, cd

import requests
import signal
from contextlib import contextmanager
from pytest import mark

slow = mark.slow


#https://www.jujens.eu/posts/en/2018/Jun/02/python-timeout-function/


@contextmanager
def timeout(time):
    # Register a function to raise a TimeoutError on the signal.
    signal.signal(signal.SIGALRM, raise_timeout)
    # Schedule the signal to be sent after ``time``
    signal.alarm(time)

    try:
        yield
    except requests.Timeout as err:
        pass
    finally:
        # Unregister the signal so it won't be triggered
        # if the timeout is not reached.
        signal.signal(signal.SIGALRM, signal.SIG_IGN)


def raise_timeout(signum, frame):
    raise TimeoutError

# def raxml_files_with_timeout():
#     # Add a timeout block.




#@slow
def test_run_raxml():

	workdir = "tests/output/test_run_raxml"
	absworkdir = os.path.abspath(workdir)
	conf = ConfigObj("tests/data/test.config", interactive=False)


	#load data 
	data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
	data_obj.workdir = absworkdir
	ids = IdDicts(conf, workdir=data_obj.workdir)
	ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

	scraper = PhyscraperScrape(data_obj, ids)
	blast_dir = "tests/data/precooked/fixed/tte_blast_files"

	# run needed functions
	# scraper.run_blast_wrapper()
	scraper.read_blast_wrapper(blast_dir=blast_dir)
	scraper.remove_identical_seqs()

	scraper.data.write_papara_files()
	scraper.align_query_seqs()
	scraper.place_query_seqs()
	scraper.est_full_tree()
	# scraper.generate_streamed_alignment()
	assert os.path.exists("{}/RAxML_bestTree.{}".format(scraper.workdir, scraper.date))
	# scraper.generate_streamed_alignment()
	if not os.path.exists("{}/previous_run".format(scraper.workdir)):
		os.mkdir("{}/previous_run".format(scraper.workdir))
	os.system("mv {}/papara_alignment.extended  {}/previous_run/papara_alignment.extended".format(scraper.workdir, scraper.workdir))

	# if os.path.exists("{}/RAxML_bootstrap.all{}".format(scraper.workdir, scraper.date)):
	# 	fn = "{}/RAxML_bootstrap.all{}".format(scraper.workdir, scraper.date)
	# 	os.system("mv {}  ./previous_run/{}".format(fn, fn))

	# CWD = os.getcwd()
	# print(CWD)
	# os.chdir(scraper.workdir)
	# with timeout(30):
	# 	print('entering block')
	# 	scraper.calculate_bootstrap()

	# 	# scraper.generate_streamed_alignment()
	# 	import time
	# 	time.sleep(50)
	# 	print('This should never get printed because the line before timed out')
	# # os.chdir(CWD)
	# # scraper.calculate_boostrap()    
	# assert os.path.exists("{}/RAxML_bootstrap.all{}".format(scraper.workdir, scraper.date))


@mark.xfail
def test_mpi():    
    env_var = [os.environ.get('PMI_RANK'), os.environ.get('PMI_SIZE'), os.environ.get('OMPI_COMM_WORLD_SIZE')]
    mpi = False
    for var in env_var:
        if var is not None:
            mpi = True
    assert mpi == True

@mark.xfail
def test_internal_mpi():
	import pickle
	import sys
	import os
	import subprocess
	from physcraper import ConfigObj, PhyscraperScrape, IdDicts
	from mpi4py import MPI

	# set up until test
	workdir = "tests/output/test_mpi_raxml"
	absworkdir = os.path.abspath(workdir)
	conf = ConfigObj("tests/data/test.config", interactive=False)


	#load data 
	data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
	data_obj.workdir = absworkdir
	ids = IdDicts(conf, workdir=data_obj.workdir)
	ids.acc_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_acc_map.p", "rb"))

	scraper = PhyscraperScrape(data_obj, ids)
	blast_dir = "tests/data/precooked/fixed/tte_blast_files"

	# run needed functions
	scraper.read_blast_wrapper(blast_dir=blast_dir)
	scraper.remove_identical_seqs()

	scraper.data.write_papara_files()
	scraper.align_query_seqs()
	scraper.place_query_seqs()
	scraper.est_full_tree()


	# scraper.generate_streamed_alignment()
	assert os.path.exists("{}/RAxML_bestTree.{}".format(scraper.workdir, scraper.date))
	# scraper.generate_streamed_alignment()
	if not os.path.exists("{}/previous_run".format(scraper.workdir)):
		os.mkdir("{}/previous_run".format(scraper.workdir))
	os.system("mv {}/papara_alignment.extended  {}/previous_run/papara_alignment.extended".format(scraper.workdir, scraper.workdir))


	cwd = os.getcwd()
	# os.chdir(scraper.workdir)


	ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
	nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
	print(nnodes, ntasks)
	env_var = int(nnodes) * int(ntasks)
	#env_var = os.environ.get('SLURM_JOB_CPUS_PER_NODE', 7)
	print(env_var)

	assert os.path.exists("{}/previous_run/papara_alignment.extended".format(scraper.workdir))
	with cd(scraper.workdir):
		print("run with mpi")
		subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxmlHPC-MPI-AVX2", 
		                 "-m", "GTRCAT",
		                 "-s", "{}/previous_run/papara_alignment.extended".format(scraper.workdir),
		                 "-p", "1", "-f", "a", "-x", "1", "-#", "autoMRE",
		                 "-n", "all{}".format(scraper.date)])
	# os.chdir(cwd)
