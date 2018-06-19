if [ -f output/owndata/ ]; then echo 'Found some!'; fi

python tests/remove_identical_seqs.py
python tests/fromfile.py
python tests/configread.py
python tests/test_species_translation.py
python tests/test_write_labelled.py
