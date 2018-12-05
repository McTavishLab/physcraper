if [ -f output/owndata/ ]; then echo 'Found some!'; fi

py.test  --cov=physcraper --runslow -v tests/test_*.py

