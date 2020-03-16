
# Get the OpenTreeeOfLife identifier from a named clade
 
from peyotl.api import APIWrapper

import sys



input = sys.argv[1]
tx = APIWrapper().taxomachine
context = None

nms = tx.TNRS([input])
for i in nms['results']:
    for j in i['matches']:
    	print j[u'taxon'][u'name']
        print j[u'taxon'][u'ott_id']