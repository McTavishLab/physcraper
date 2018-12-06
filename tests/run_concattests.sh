if [ -f output/owndata/ ]; then echo 'Found some!'; fi

py.test --runslow --cov=physcraper tests/test_concat_make_concat_id.py
py.test --runslow --cov=physcraper tests/test_concat_select_rnd_seq.py
py.test --runslow --cov=physcraper tests/test_concat_spcounter.py
py.test --runslow --cov=physcraper tests/test_concat_sp_to_keep.py
py.test --runslow --cov=physcraper tests/test_concat_run.py
