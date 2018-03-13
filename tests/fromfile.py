from physcraper import generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape


#Use OpenTree phylesystem identifiers to get study and tree
seqaln = "tests/data/input.fas"
mattype="fasta"
workdir="tests/fromfile"
configfi = "example.config"
treefile = "tests/data/input.tre"
otu_json = "tests/data/otu_dict.json"

data_obj = generate_ATT_from_files(seqaln,
                        mattype,
                        workdir,
                        treefile,
                        otu_json)


conf = ConfigObj(configfi)
ids = IdDicts(conf, workdir=workdir)
scraper = PhyscraperScrape(data_obj, ids, conf)
