# package import
from physcraper import get_mrca_ott



l1 = ['515698','590452','643717']
taxon_not_in_tree = ['372706','563165']

def test_get_mrca_ott():
    a1 = get_mrca_ott(l1)
    a2 = get_mrca_ott(taxon_not_in_tree)
    assert a1 == 1042120


