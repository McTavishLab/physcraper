import subprocess

def test_physcraper_run_script():
    subprocess.check_call(["python", "bin/physcraper_run.py", 
                            "-s", "ot_350", 
                            "-t", "Tr53297",
                            "-a", "docs/examples/ot_350Tr53297.aln",
                            "-as", "nexus",
                            "-no_est",
                            "-o", "tests/data/precooked/ot_350/"])
