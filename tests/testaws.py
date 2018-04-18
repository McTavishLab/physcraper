import physcraper.AWSWWW as AWSWWW
import sys
import datetime
from physcraper import ConfigObj

query = 'TTGACCTCGGATCAGGTAGGAATACCCGCTGAACTTAAGCATATCAATAAGCGGAGGAAAAGAAACCAACAGGGATTGCCCCAGTAACGGCGAGTGAAGCGGCAACAGCTCAAATTTGAAATCTGGCCCCAGGCCCGAGTTGTAATTTGCAGAGGATGCTTCGGGCGCGACGCCTTCCAAGTCCCCTGGAACGGGGCGCCTTAGAGGGTGAGAGCCCCGTACGGTTGGACGTCAAGCCTGTGTGAAGCTCCTTCGACGAGTCGAGTAGTTTGGGAATGCTGCTCAAAATGGGAGGTAGACCCCTTCTAAAGCTAAATACCGGCCAGAGACCGATAGCGCACAAGTAGAGTGATCGAAAGATGAAAAGCACTTTGAAAAGAGGGTTAAACAGCACGTGAAATTGTTGAAAGGGAAGCGCTCGTGACCAGACTTGCGCCGGGGCGATCATCTGGCGTTCTCGCCGGTGCACTCGCCCCGGCTCAGGCCAGCGTCGGTTCGGGAGGGGGGACAAAGGCGTCGGGGATGTGGCTCCCTCGGGAGTGTTATAGCCCGGCGTGCAATGCCCCCGCCCGGACCGAGGTTCGCGCTCTGCAAGGACGCTGGCGTAATGGTCACCAGCGGCCCGTCTTGAAACACGGACCAAGGAGTCGAGGTTCTGCGCGAGTGTTTGGGTGTCAAACCCGCACGCGTAATGAAAGTGAACGTAGGTGGGAGCCTCGGCGCACCACCGACCGATCCTGATGTTCTCGGATGGATTTGAGTAGGAGCGTAGGGCCTCGGACCCGAAAGATGGTGAACTATGCTTGGATAGGGTGAAGCCAGAGGAAACTCTGGTGGAGGCTCGCAGCGGTTCTGACGTGCAAATCGATCGTCAAATCTGAGCATGGGGGCGAAAGACTAATCGAACCATCTAGTAGCTGGTTACCGCCGAAGTTTCCCTCAGGATAGCAGTGTT'


configfi = "tests/data/aws.config"

conf = ConfigObj(configfi)


start = datetime.datetime.now()
sys.stderr.write("running a blast query against AWS instance START {}".format(start))



result_handle = AWSWWW.qblast("blastn",
                               "nt",
                               query,
                               url_base = conf.url_base,
                               hitlist_size=2,
                               num_threads=4)


end = datetime.datetime.now()
sys.stderr.write("END {}".format(end))
sys.stderr.write(result_handle.read())
result_handle.close()

start = datetime.datetime.now()
sys.stderr.write("running a blast query against NCBI START {}".format(start))



result_handle = AWSWWW.qblast("blastn",
                               "nt",
                               query,
                               hitlist_size=2)


end = datetime.datetime.now()
sys.stderr.write("END {}".format(end))
sys.stderr.write(result_handle.read())
result_handle.close()