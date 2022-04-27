
import physcraper
from physcraper import viz
from physcraper import TreeTax

t_file = "pg_55_test/outputs_pg_55tree5864/labelled_pg_55tree5864.tre"
csv_file = "pg_55_test/outputs_pg_55tree5864/otu_info_pg_55tree5864.csv"


physcraper.viz.plot_tree(tree_file = t_file,
                         pdf_file = "pg_55_test/labelled_pg_55tree5864.pdf",
                         font_size = 8)


