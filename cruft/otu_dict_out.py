import pickle
otu_dict = pickle.load(open('otu_dict.p','rb'))

fi=open('tip_info.txt','w')


info_cats = ['^ot:ottTaxonName', '^ot:treebaseOTUId', '^physcraper:status', u'^ot:ottId', 'physcraper:status', '^ot:originalLabel', '^ncbi:taxon', '^ncbi:title', '^ncbi:accession', '^ncbi:gi', '^physcraper:last_blasted']

fi.write("\t".join(['uniqueID']+info_cats))
fi.write("\n")

for tip in otu_dict:
    fi.write("\t".join([tip]+[str(otu_dict[tip].get(cat)) for cat in info_cats]))
    fi.write("\n")

fi.close()