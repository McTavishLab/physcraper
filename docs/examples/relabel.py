
import physcraper
from physcraper import TreeTax

otu_t_file = "pg_55_web/outputs_pg_55tree5864/physcraper_pg_55tree5864.tre"
otu_info_file = "pg_55_web/run_pg_55tree5864/otu_info_pg_55tree5864.json"


tt = TreeTax(otu_json=otu_info_file, 
             treefrom=otu_t_file )



tt.write_labelled(path="ids_only.tre",
                  label = "^ot:ottId",
                  norepeats=True)



