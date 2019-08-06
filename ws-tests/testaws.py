import physcraper.AWSWWW as AWSWWW
from Bio.Blast import NCBIXML
import sys
import datetime
from physcraper import ConfigObj

query = 'TTGACCTCGGATCAGGTAGGAATACCCGCTGAACTTAAGCATATCAATAAGCGGAGGAAAAGAAACCAACAGGGATTGCCCCAGTAACGGCGAGTGAAGCGGCAACAGCTCAAATTTGAAATCTGGCCCCAGGCCCGAGTTGTAATTTGCAGAGGATGCTTCGGGCGCGACGCCTTCCAAGTCCCCTGGAACGGGGCGCCTTAGAGGGTGAGAGCCCCGTACGGTTGGACGTCAAGCCTGTGTGAAGCTCCTTCGACGAGTCGAGTAGTTTGGGAATGCTGCTCAAAATGGGAGGTAGACCCCTTCTAAAGCTAAATACCGGCCAGAGACCGATAGCGCACAAGTAGAGTGATCGAAAGATGAAAAGCACTTTGAAAAGAGGGTTAAACAGCACGTGAAATTGTTGAAAGGGAAGCGCTCGTGACCAGACTTGCGCCGGGGCGATCATCTGGCGTTCTCGCCGGTGCACTCGCCCCGGCTCAGGCCAGCGTCGGTTCGGGAGGGGGGACAAAGGCGTCGGGGATGTGGCTCCCTCGGGAGTGTTATAGCCCGGCGTGCAATGCCCCCGCCCGGACCGAGGTTCGCGCTCTGCAAGGACGCTGGCGTAATGGTCACCAGCGGCCCGTCTTGAAACACGGACCAAGGAGTCGAGGTTCTGCGCGAGTGTTTGGGTGTCAAACCCGCACGCGTAATGAAAGTGAACGTAGGTGGGAGCCTCGGCGCACCACCGACCGATCCTGATGTTCTCGGATGGATTTGAGTAGGAGCGTAGGGCCTCGGACCCGAAAGATGGTGAACTATGCTTGGATAGGGTGAAGCCAGAGGAAACTCTGGTGGAGGCTCGCAGCGGTTCTGACGTGCAAATCGATCGTCAAATCTGAGCATGGGGGCGAAAGACTAATCGAACCATCTAGTAGCTGGTTACCGCCGAAGTTTCCCTCAGGATAGCAGTGTT'


configfi = "tests/data/aws.config"

#try:
conf = ConfigObj(configfi)

start = datetime.datetime.now()
sys.stderr.write("running a blast query against AWS instance START {}".format(start))

#conf.url_base = "https://ec2-54-175-193-21.compute-1.amazonaws.com/Blast.cgi"


conf.url_base = "https://ec2-54-175-193-21.compute-1.amazonaws.com/cgi-bin/blast.cgi"
result_handle = AWSWWW.qblast("blastn",
                                 "nt",
                                 query,
                                 url_base = conf.url_base,
                                 hitlist_size=2,
                                 num_threads=4)

result_handle = AWSWWW.qblast("blastn",
                                 "nt",
                                 query,
                                 hitlist_size=2,
                                 num_threads=4)


end = datetime.datetime.now()
sys.stderr.write("END {}".format(end))


blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
   assert len(blast_record.alignments) == 2
                                   
result_handle.close()

