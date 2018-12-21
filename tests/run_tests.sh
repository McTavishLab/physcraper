if [ -f output/owndata/ ]; then echo 'Found some!'; fi

py.test  --cov=physcraper --runslow  --runconcat -v tests/


#py.test --cov=physcraper tests/test_remove_identical_seqs.py
# py.test  --cov=physcraper --runslow -v \
#  tests/test_remove_identical_seqs.py \
#  tests/test_fromfile.py \
#  tests/test_configread.py \
#  tests/test_species_translation.py \
#  tests/test_write_labelled.py \
#  tests/test_add_all.py \
#  tests/test_loop_for_blast.py \
#  tests/test_prune_short.py \
#  tests/test_remove_taxa_aln_tre.py \
#  tests/test_run_local_blast.py \
#  tests/test_sp_d.py \
#  tests/test_sp_seq_d.py \
#  tests/test_select_seq_by_local_blast.py \
#  tests/test_calc_mean_sd.py \
#  tests/test_read_local_blast.py \
#  tests/tests_write_blast.py \
#  tests/test_remove_id_seq.py \
#  tests/test_addLocal.py \
#  tests/test_reconcile.py \
#  tests/test_trim.py \
#  tests/test_unmapped_taxa.py \
#  tests/test_blacklist.py
