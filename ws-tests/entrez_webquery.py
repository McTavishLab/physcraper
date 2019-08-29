# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 09:49:13 2018

@author: Martha Kandziora
"""


### Error message
#  File "/home/mkandziora/physcrape_web_dev/physcraper/physcraper/wrappers.py", line 560, in filter_data_run
#    filteredScrape.run_blast(delay=14)
#  File "/home/mkandziora/physcrape_web_dev/physcraper/physcraper/__init__.py", line 1930, in run_blast
#    self.run_web_blast_query(query, equery, fn_path)
#  File "/home/mkandziora/physcrape_web_dev/physcraper/physcraper/__init__.py", line 1831, in run_web_blast_query
#    hitlist_size=self.config.hitlist_size)
#  File "/home/mkandziora/physcrape_web_dev/physcraper/physcraper/AWSWWW.py", line 142, in qblast
#    handle = _urlopen(request)
#  File "/usr/lib64/python2.7/urllib2.py", line 154, in urlopen
#    return opener.open(url, data, timeout)
#  File "/usr/lib64/python2.7/urllib2.py", line 431, in open
#    response = self._open(req, data)
#  File "/usr/lib64/python2.7/urllib2.py", line 449, in _open
#    '_open', req)
#  File "/usr/lib64/python2.7/urllib2.py", line 409, in _call_chain
#    result = func(*args)
#  File "/usr/lib64/python2.7/urllib2.py", line 1258, in https_open
#    context=self._context, check_hostname=self._check_hostname)
#  File "/usr/lib64/python2.7/urllib2.py", line 1214, in do_open
#    raise URLError(err)
#urllib2.URLError: <urlopen error [Errno -2] Name or service not known>

## from printing statements
#use BLAST webservice
#TCGAAACCTGCATAGCAGAACGACCCGTGAACATGTAACAACAACCGGGTGTCCTTGGTATCGGACTCTTGTCTGATTCTTTGGATGCCTCGTCGATGTGCGTCTTTGGCCTGCCCCTTGGGTTCCAATTACGTCACGTTGGCACAACAACAACCCCCCGGCACGGCATGTGCCAAGGAAATTTAAACATAAGAAGGGCTCGTACCATGCATCCCCGTTCGCGGGGTTTGCGTGGGATGTGGCTTCTTTATAATCACAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTCGAGGGCACGTCTGTCTGGGCGTCACACATCGCGTCGCCCCCACCATTCCTCTCTGATCGGGATGATGGGAAGGGGGCGGATATTGGTCTCCCGTTCCTACGGTGCGGTTGGCTAAAATAGGAGTCCCCTTCGACGGACGCACGACTAGTGGTGGTTGACAAGACCCTCTTATCGAGTCGTGCGTTCTAAGGAGTAAGGAAGATCTCTTAAAAAACCC
#txid102812[orgn] AND 1900/01/01:2018/10/23[mdat]
#5000

import physcraper.AWSWWW as AWSWWW


query = "TCGAAACCTGCATAGCAGAACGACCCGTGAACATGTAACAACAACCGGGTGTCCTTGGTATCGGACTCTTGTCTGATTCTTTGGATGCCTCGTCGATGTGCGTCTTTGGCCTGCCCCTTGGGTTCCAATTACGTCACGTTGGCACAACAACAACCCCCCGGCACGGCATGTGCCAAGGAAATTTAAACATAAGAAGGGCTCGTACCATGCATCCCCGTTCGCGGGGTTTGCGTGGGATGTGGCTTCTTTATAATCACAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGTCGAGGGCACGTCTGTCTGGGCGTCACACATCGCGTCGCCCCCACCATTCCTCTCTGATCGGGATGATGGGAAGGGGGCGGATATTGGTCTCCCGTTCCTACGGTGCGGTTGGCTAAAATAGGAGTCCCCTTCGACGGACGCACGACTAGTGGTGGTTGACAAGACCCTCTTATCGAGTCGTGCGTTCTAAGGAGTAAGGAAGATCTCTTAAAAAACCC"
equery = "txid102812[orgn]"

result_handle = AWSWWW.qblast("blastn",
                              "nt",
                              query,
                              entrez_query=equery,
                              hitlist_size=50)


tmp = result_handle.read()
result_handle.close()


