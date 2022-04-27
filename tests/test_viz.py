
import physcraper
from physcraper import viz

t_file = "pg_55_test/outputs_pg_55tree5864/labelled_pg_55tree5864.tre"

physcraper.viz.plot_tree(tree_file = t_file,
                    pdf_file = "pg_55_test/labelled_pg_55tree5864.pdf",
                    font_size = 8)
