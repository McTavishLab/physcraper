import sys
import os
import json
import random
import pickle
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
from physcraper import ConfigObj, debug, IdDicts, FilterBlast
#


#
seqaln= "tests/data/Senecio_its_input/senecio_its.fasta"
mattype="fasta"
trfn= "tests/data/Senecio_its_input/its_new.tre"
schema_trf = "newick"
id_to_spn = r"tests/data/Senecio_its_input/uniquetip_to_name_its.csv"

workdir="Senecio__alientaxa"
configfi = "tests/data/blubb_localblast.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold = 2
selectby = "blast"
downtorank = None
# downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None
blast_dir = "{}/save_xml".format(workdir)






random.seed(1234)
try:
	if os.path.isfile("{}/Senecio_alien_round2_crashpoint.p".format(workdir)): 
		filteredScrape = pickle.load(open("{}/Senecio_alien_round2_crashpoint.p".format(workdir),'rb'))
		filteredScrape.generate_streamed_alignment(treshold)
		filteredScrape.dump("{}/Senecio_alien_round2.p".format(workdir))

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)
	else:

	# round 3
		if os.path.isfile("{}/Senecio_alien_round3.p".format(workdir)): 
			sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
			filteredScrape = pickle.load(open("{}/Senecio_alien_round3.p".format(workdir),'rb'))
			filteredScrape.repeat = 1   
			filteredScrape.read_blast(blast_dir=blast_dir)

		elif filteredScrape.repeat == 1: 
			filteredScrape.run_blast()
			filteredScrape.read_blast(blast_dir=blast_dir)
			filteredScrape.remove_identical_seqs()

			folder = '{}/blast/'.format(filteredScrape.workdir)
			for the_file in os.listdir(folder):
				file_path = os.path.join(folder, the_file)
				if os.path.isfile(file_path):
					os.unlink(file_path)
			debug("make sp_dict")	
			if treshold != None:  
				filteredScrape.sp_dict(downtorank)
				filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
				filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
				filteredScrape.replace_new_seq()
			filteredScrape.generate_streamed_alignment(treshold)
			filteredScrape.dump("{}/Senecio_alien_round3.p".format(workdir))
			#debug_count+=1
			#filteredScrape.dump('{}/round{}.p'.format(workdir, debug_count))

			folder = '{}/blast/'.format(filteredScrape.workdir)
			for the_file in os.listdir(folder):
				file_path = os.path.join(folder, the_file)
				if os.path.isfile(file_path):
					os.unlink(file_path)


	# round 4
		if os.path.isfile("{}/Senecio_alien_round4.p".format(workdir)): 
			sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
			filteredScrape = pickle.load(open("{}/Senecio_alien_round4.p".format(workdir),'rb'))
			filteredScrape.repeat = 1   

		elif filteredScrape.repeat == 1: 
			filteredScrape.run_blast()
			filteredScrape.read_blast(blast_dir=blast_dir)
			filteredScrape.remove_identical_seqs()

			folder = '{}/blast/'.format(filteredScrape.workdir)
			for the_file in os.listdir(folder):
				file_path = os.path.join(folder, the_file)
				if os.path.isfile(file_path):
					os.unlink(file_path)
			debug("make sp_dict")	
			if treshold != None:  
				filteredScrape.sp_dict(downtorank)
				filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
				filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
				filteredScrape.replace_new_seq()
			filteredScrape.generate_streamed_alignment(treshold)
			filteredScrape.dump("{}/Senecio_alien_round4.p".format(workdir))
			#debug_count+=1
			#filteredScrape.dump('{}/round{}.p'.format(workdir, debug_count))

			folder = '{}/blast/'.format(filteredScrape.workdir)
			for the_file in os.listdir(folder):
				file_path = os.path.join(folder, the_file)
				if os.path.isfile(file_path):
					os.unlink(file_path)


#################################################################
except: # look for other dumped files....			

	if os.path.exists(otu_jsonfi):
		otu_json = json.load(open(otu_jsonfi))
	else:
		otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
		if not os.path.exists(workdir):
		   os.mkdir(workdir)
		json.dump(otu_json, open(otu_jsonfi,"w"))


	if os.path.isfile("{}/Senecio_alien_round1.p".format(workdir)): 
		sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
		filteredScrape = pickle.load(open("{}/Senecio_alien_round1.p".format(workdir),'rb'))
		filteredScrape.repeat = 1   

	else:  

		sys.stdout.write("setting up Data Object\n")
		sys.stdout.flush()
		#read the config file into a configuration object
		conf = ConfigObj(configfi)
		# print("config")
		debug(dir(conf))
		debug(conf.email)


		#Generate an linked Alignment-Tree-Taxa object
		data_obj = generate_ATT_from_files(seqaln=seqaln, 
							 mattype=mattype, 
							 workdir=workdir,
							 treefile=trfn,
							 schema_trf=schema_trf,
							 otu_json=otu_jsonfi,
							 #email=conf.email,
							 ingroup_mrca=None)

		#Prune sequnces below a certain length threshold
		#This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
		data_obj.prune_short()
		data_obj.write_files()

		data_obj.write_labelled( label='user:TaxonName')
		data_obj.write_otus("otu_info", schema='table')
		data_obj.dump()

		#ids = IdDicts(conf, workdir="example")
		#		 if os.path.isfile("{}/id_pickle.p".format(workdir)): 

		#		 #if os.path.isfile(conf.id_pickle):
		#			 sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
		# #		thawed_id = open(conf.id_json, 'r').readlines()
		# #		ids = jsonpickle.decode(thawed_id)
		# #		scraper.repeat = 1
		#			 ids = pickle.load(open("{}/id_pickle.p".format(workdir),'rb'))
		#		 else:
		sys.stdout.write("setting up id dictionaries\n")
		sys.stdout.flush()
		# if os.path.isfile("{}/id_pickle.p".format(workdir)): 
		#	 sys.stdout.write("Reloading from pickled scrapefile: id\n")
		#	 ids = pickle.load(open("{}/id_pickle.p".format(workdir),'rb'))

		# else:   
		ids = IdDicts(conf, workdir=workdir)
		ids.dump()

		# if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
		#	 sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
		#	 scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
		#	 scraper.repeat = 1	
		# else:   
			#Now combine the data, the ids, and the configuration into a single physcraper scrape object
		filteredScrape =  FilterBlast(data_obj, ids)
		if add_local_seq != None:
			debug("will add local sequences now")
			filteredScrape.add_local_seq(add_local_seq, id_to_spn_addseq_json)
			# scraper.replace_new_seq()
			filteredScrape.remove_identical_seqs()
			filteredScrape.generate_streamed_alignment(treshold)
		#run the ananlyses
		filteredScrape.run_blast()
		filteredScrape.read_blast(blast_dir=blast_dir)
		filteredScrape.remove_identical_seqs()
		filteredScrape.dump()
		debug(treshold)
		if treshold != None:  
			filteredScrape.sp_dict(downtorank)
			filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
			filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
			filteredScrape.replace_new_seq()
		debug("from replace to streamed aln")
		filteredScrape.generate_streamed_alignment(treshold)
		filteredScrape.dump("{}/Senecio_alien_round1.p".format(workdir))
		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)

	# round 2
	if os.path.isfile("{}/Senecio_alien_round2.p".format(workdir)): 
		sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
		filteredScrape = pickle.load(open("{}/Senecio_alien_round2.p".format(workdir),'rb'))
		filteredScrape.repeat = 1   
	elif filteredScrape.repeat == 1: 
		filteredScrape.run_blast()
		filteredScrape.read_blast(blast_dir=blast_dir)
		filteredScrape.remove_identical_seqs()

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)
		debug("make sp_dict")	
		if treshold != None:  
			filteredScrape.sp_dict(downtorank)
			filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
			filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
			filteredScrape.replace_new_seq()
			filteredScrape.dump("{}/Senecio_alien_round2_crashpoint.p".format(workdir))

		filteredScrape.generate_streamed_alignment(treshold)
		filteredScrape.dump("{}/Senecio_alien_round2.p".format(workdir))
		#debug_count+=1
		#filteredScrape.dump('{}/round{}.p'.format(workdir, debug_count))

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)

	# round 3
	if os.path.isfile("{}/Senecio_alien_round3.p".format(workdir)): 
		sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
		filteredScrape = pickle.load(open("{}/Senecio_alien_round3.p".format(workdir),'rb'))
		filteredScrape.repeat = 1   

	elif filteredScrape.repeat == 1: 
		filteredScrape.run_blast()
		filteredScrape.read_blast(blast_dir=blast_dir)
		filteredScrape.remove_identical_seqs()

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)
		debug("make sp_dict")	
		if treshold != None:  
			filteredScrape.sp_dict(downtorank)
			filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
			filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
			filteredScrape.replace_new_seq()
		filteredScrape.generate_streamed_alignment(treshold)
		filteredScrape.dump("{}/Senecio_alien_round3.p".format(workdir))
		#debug_count+=1
		#filteredScrape.dump('{}/round{}.p'.format(workdir, debug_count))

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)


	# round 4
	if os.path.isfile("{}/Senecio_alien_round4.p".format(workdir)): 
		sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
		filteredScrape = pickle.load(open("{}/Senecio_alien_round4.p".format(workdir),'rb'))
		filteredScrape.repeat = 1   

	elif filteredScrape.repeat == 1: 
		filteredScrape.run_blast()
		filteredScrape.read_blast(blast_dir=blast_dir)
		filteredScrape.remove_identical_seqs()

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)
		debug("make sp_dict")	
		if treshold != None:  
			filteredScrape.sp_dict(downtorank)
			filteredScrape.make_sp_seq_dict(treshold=treshold, selectby=selectby)
			filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
			filteredScrape.replace_new_seq()
		filteredScrape.generate_streamed_alignment(treshold)
		filteredScrape.dump("{}/Senecio_alien_round4.p".format(workdir))
		#debug_count+=1
		#filteredScrape.dump('{}/round{}.p'.format(workdir, debug_count))

		folder = '{}/blast/'.format(filteredScrape.workdir)
		for the_file in os.listdir(folder):
			file_path = os.path.join(folder, the_file)
			if os.path.isfile(file_path):
				os.unlink(file_path)