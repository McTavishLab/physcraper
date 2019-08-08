import physcraper.AWSWWW as AWSWWW
from Bio.Blast import NCBIXML
import sys
import time
from physcraper import ConfigObj

query = 'TTGACCTCGGATCAGGTAGGAATACCCGCTGAACTTAAGCATATCAATAAGCGGAGGAAAAGAAACCAACAGGGATTGCCCCAGTAACGGCGAGTGAAGCGGCAACAGCTCAAATTTGAAATCTGGCCCCAGGCCCGAGTTGTAATTTGCAGAGGATGCTTCGGGCGCGACGCCTTCCAAGTCCCCTGGAACGGGGCGCCTTAGAGGGTGAGAGCCCCGTACGGTTGGACGTCAAGCCTGTGTGAAGCTCCTTCGACGAGTCGAGTAGTTTGGGAATGCTGCTCAAAATGGGAGGTAGACCCCTTCTAAAGCTAAATACCGGCCAGAGACCGATAGCGCACAAGTAGAGTGATCGAAAGATGAAAAGCACTTTGAAAAGAGGGTTAAACAGCACGTGAAATTGTTGAAAGGGAAGCGCTCGTGACCAGACTTGCGCCGGGGCGATCATCTGGCGTTCTCGCCGGTGCACTCGCCCCGGCTCAGGCCAGCGTCGGTTCGGGAGGGGGGACAAAGGCGTCGGGGATGTGGCTCCCTCGGGAGTGTTATAGCCCGGCGTGCAATGCCCCCGCCCGGACCGAGGTTCGCGCTCTGCAAGGACGCTGGCGTAATGGTCACCAGCGGCCCGTCTTGAAACACGGACCAAGGAGTCGAGGTTCTGCGCGAGTGTTTGGGTGTCAAACCCGCACGCGTAATGAAAGTGAACGTAGGTGGGAGCCTCGGCGCACCACCGACCGATCCTGATGTTCTCGGATGGATTTGAGTAGGAGCGTAGGGCCTCGGACCCGAAAGATGGTGAACTATGCTTGGATAGGGTGAAGCCAGAGGAAACTCTGGTGGAGGCTCGCAGCGGTTCTGACGTGCAAATCGATCGTCAAATCTGAGCATGGGGGCGAAAGACTAATCGAACCATCTAGTAGCTGGTTACCGCCGAAGTTTCCCTCAGGATAGCAGTGTT'


configfi = "tests/data/aws.config"

#try:
conf = ConfigObj(configfi)



#conf.url_base = "https://ec2-54-175-193-21.compute-1.amazonaws.com/Blast.cgi"

threads = 20


start = time.time()
result_handle = AWSWWW.qblast("blastn",
                                 "nt",
                                 query,
                                 url_base = conf.url_base,
                                 hitlist_size=2,
                                 num_threads=threads)

end = time.time()
sys.stderr.write("AWS: Threads {} Time  {}\n".format(threads, end-start))

start = time.time()
sys.stderr.write("running a blast query against NCBI  START {}\n".format(start))
result_handle = AWSWWW.qblast("blastn",
                                     "nt",
                                     query,
                                     hitlist_size=2)
end = time.time()
sys.stderr.write("NCBI: Time  {}\n".format(end-start))




blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
   assert len(blast_record.alignments) == 2
                                   
result_handle.close()

