import physcraper
from physcraper import TreeTax
from physcraper import viz




otu_t_file = "pg_55_web/outputs_pg_55tree5864/physcraper_pg_55tree5864.tre"
otu_info_file = "pg_55_web/run_pg_55tree5864/otu_info_pg_55tree5864.json"


tt = TreeTax(otu_json=otu_info_file, 
             treefrom=otu_t_file )

new_tip_labels = []
for tip in tt.otu_dict:
    taxon_name = tt.otu_dict[tip]['^ot:ottTaxonName']
    if tt.otu_dict[tip]['^physcraper:status'] == 'new':
        new_tip_labels.append(taxon_name.replace(' ', '_')+'_'+tip)


print(new_tip_labels)

t_file = "pg_55_web/outputs_pg_55tree5864/labelled_pg_55tree5864.tre"

physcraper.viz.plot_tree(tree_file = t_file,
                         pdf_file = "testing_labelling.pdf",
                         font_size = 8,
                         list_of_new_taxa = new_tip_labels)


