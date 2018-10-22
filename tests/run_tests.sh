if [ -f output/owndata/ ]; then echo 'Found some!'; fi

python tests/test_remove_identical_seqs.py
python tests/test_fromfile.py
python tests/test_configread.py
python tests/test_reconcile.py
python tests/test_species_translation.py
python tests/test_write_labelled.py
python tests/test_add_all.py
python tests/test_loop_for_blast.py
python tests/test_prune_short.py
python tests/test_remove_taxa_aln_tre.py
python tests/test_run_local_blast.py
python tests/test_sp_d.py
python tests/test_sp_seq_d.py
python tests/test_select_seq_by_local_blast.py
python tests/test_calc_mean_sd.py
python tests/test_read_local_blast.py
python tests/tests_write_blast.py
python tests/test_remove_id_seq.py
python tests/test_addLocal.py
python tests/test_blacklist.py
