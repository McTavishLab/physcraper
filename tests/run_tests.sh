if [ -f output/owndata/ ]; then echo 'Found some!'; fi

py.test tests/test_remove_identical_seqs.py
py.test tests/test_fromfile.py
py.test tests/test_configread.py
py.test tests/test_species_translation.py
py.test tests/test_write_labelled.py
py.test tests/test_add_all.py
py.test tests/test_loop_for_blast.py
py.test tests/test_prune_short.py
py.test tests/test_remove_taxa_aln_tre.py
py.test tests/test_run_local_blast.py
py.test tests/test_sp_d.py
py.test tests/test_sp_seq_d.py
py.test tests/test_select_seq_by_local_blast.py
py.test tests/test_calc_mean_sd.py
py.test tests/test_read_local_blast.py
py.test tests/tests_write_blast.py
py.test tests/test_remove_id_seq.py
py.test tests/test_addLocal.py
py.test tests/test_reconcile.py
py.test tests/test_trim.py
py.test tests/test_unmapped_taxa.py
python tests/test_trim2.py

py.test tests/test_blacklist.py
