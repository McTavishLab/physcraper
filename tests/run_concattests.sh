if [ -f output/owndata/ ]; then echo 'Found some!'; fi

python tests/test_concat_make_concat_id.py
python tests/test_concat_select_rnd_seq.py
python tests/test_concat_spcounter.py
python tests/test_concat_sp_to_keep.py
