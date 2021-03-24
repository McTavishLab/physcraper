import subprocess

def test_physcraper_run_script():
    subprocess.check_call(["python", "bin/physcraper_run.py", 
                            "-s", "ot_350", 
                            "-t", "Tr53297",
                            "-a", "docs/examples/inputdata/ot_350Tr53297.aln",
                            "-as", "nexus",
                            "-no_est",
                            "-o", "tests/data/precooked/ot_350/"])



def test_physcraper_rerun_script():
    subprocess.check_call(["python", "bin/physcraper_run.py", 
                            "-re", "tests/data/precooked/ot_350",
                            "-o", "tests/tmp/ot350",
                            "-tag", "rerun",
                            "-no_est"])

def test_script_args():
    subprocess.check_call(["python", "bin/physcraper_run.py",
                            "-tf", "tests/data/tiny_test_example/test.tre",
                            "-tfs", "newick",
                            "-a", "tests/data/tiny_test_example/test.fas",
                            "--taxon_info", "tests/data/tiny_test_example/main.json",
                            "-as", "fasta",
                            "-o", "tests/tmp/owndata",
                            "-c", "tests/data/test.config",
                            "-st", "ott:1084160",
                            "-spn", "3",
                            "-bl", 'JX895264.1',
                            "-tp", "0.9",
                            "-rlmin", ".6",
                            "-rlmax", "2",
                            "-ev", '0.0001',
                            "-hl", '20',
                            "-nt", '6',
                            "-de", '30',
                            "-bs", "50",
                            "-v",
                            "-no_est"])




def test_comparison():
    subprocess.check_call(["python", "bin/tree_comparison.py",
                           "-d", "docs/examples/pg_55_local/",
                            "-og", "otu376420", "otu376439", "otu376452",
                            "-o", "tests/tmp/pg_55_comparison"])

def test_concat():
    subprocess.check_call(["python", "bbin/multi_loci.py",
                            "-d", "tests/data/precooked/multi_loc/",
                            "-f", "astral",
                            "-o", "tests/tmp/mini_astral"
                            ])
