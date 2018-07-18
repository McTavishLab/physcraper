if [ -f output/owndata/ ]; then echo 'Found some!'; fi

coverage run tests/test_remove_identical_seqs.py
coverage run tests/test_fromfile.py
coverage run tests/test_configread.py
coverage run tests/test_species_translation.py
coverage run tests/test_write_labelled.py
coverage run tests/test_add_all.py
coverage run tests/test_loop_for_blast.py
coverage run tests/test_prune_short.py
coverage run tests/test_remove_taxa_aln_tre.py
coverage run tests/test_run_local_blast.py
coverage run tests/test_sp_d.py
coverage run tests/test_sp_seq_d.py
coverage run tests/test_select_seq_by_local_blast.py
coverage run tests/test_calc_mean_sd.py
coverage run tests/test_read_local_blast.py
coverage run tests/tests_write_blast.py
coverage run tests/test_remove_id_seq.py

coverage report -m  physcraper/__init__.py physcraper/AWSWWW.py


# python tests/test_blacklist.py
