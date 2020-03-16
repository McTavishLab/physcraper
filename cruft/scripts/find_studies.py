# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 16:33:10 2018

@author: Martha Kandziora

find studies and trees from OToL
"""

from peyotl.api import APIWrapper
import sys


input_name = sys.argv[1]

a = APIWrapper()
oti = a.oti
tre = oti.find_trees(ottTaxonName=input_name)

if tre != []:   
    for item in tre:
        print("ot:studyId")
        print(item[u'ot:studyId'])
        #print(item[u'matched_trees'])
        for study in item[u'matched_trees']:
            print("tree_id")
            print(study[u'nexson_id'])

else:
    print("No study found")
