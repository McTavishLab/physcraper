
import sys
from physcraper import ConfigObj


sys.stdout.write("\nTesting configuration object contents'\n")
try:
   configfi = "tests/data/localblast.config"
   conf = ConfigObj(configfi)
   assert conf.email == todo
   assert conf.base_url == todo
   assert conf.keys == expected_keys
   sys.stdout.write("\nTest configread.py passed\n")
except:
	sys.stdout.write("\nTest configread.py FAILED (expected, test in progress)\n")