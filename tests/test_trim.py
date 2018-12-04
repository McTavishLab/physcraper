from dendropy import Tree, \
      DnaCharacterMatrix, \
      DataSet, \
      datamodel
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax, OtuJsonDict, IdDicts, ConfigObj
import os
import json

def test_trim():
  #------------------------
  seqaln= "tests/data/tiny_test_example/test_extralongseq.fas"
  mattype="fasta"
  treefile= "tests/data/tiny_test_example/test.tre"
  schema_trf = "newick"
  workdir="tests/output/test_trim"
  configfi = "tests/data/test.config"
  id_to_spn = r"tests/data/tiny_test_example/test_nicespl.csv"
  otu_jsonfi = "{}/otu_dict.json".format(workdir)



  if not os.path.exists("{}".format(workdir)):
          os.makedirs("{}".format(workdir))

  conf = ConfigObj(configfi, interactive=False)
  ids = IdDicts(conf, workdir=workdir)

  if os.path.exists(otu_jsonfi):
      print("load json")
      otu_json = json.load(open(otu_jsonfi))
  else:
      otu_json = OtuJsonDict(id_to_spn, ids)
      json.dump(otu_json, open(otu_jsonfi,"w"))


  data_obj = generate_ATT_from_files(seqaln=seqaln, 
                                   mattype=mattype, 
                                   workdir=workdir,
                                   treefile=treefile,
                                   schema_trf = schema_trf,
                                   otu_json=otu_jsonfi,
                                   ingroup_mrca=None)

  for tax, seq in data_obj.aln.items():
  	len_start = len(seq)
  	next
  data_obj.trim()
  for tax, seq in data_obj.aln.items():
  	len_end = len(seq)

  assert len_start != len_end
