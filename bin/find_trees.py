#!/usr/bin/env python
import sys
import argparse
from opentree import OT



parser = argparse.ArgumentParser()
parser.add_argument("-t","--taxon_name", help="Name of search taxon")
parser.add_argument("-ott","--ott_id", help="Name of search taxon")
parser.add_argument("-tb","--treebase", action = 'store_true', help="Rturn studies with treebase data only")
parser.add_argument("-o","--output", help="Output file path")
args = parser.parse_args()


args = parser.parse_args()
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)


assert(args.taxon_name or args.ottid), "A taxon name or an OTT id are required for search."

if args.taxon_name:
    try:
        ottid = OT.get_ottid_from_name(args.taxon_name)
    except:
        sys.stdout.write("no match to taxon name. Try finding your taxon on tree.opentreeoflife.org and inputting the taxon id using -ott\n")
        sys.exit()


if args.ott_id:
    ottid = args.ott_id

sys.stdout.write("OTT id {}\n".format(ottid))
phylesystem_studies_resp = OT.find_trees(ottid, search_property ='ot:ottId')



cites_phyl = "Members of {} present in the following studies in the OpenTree Phylesystem\n".format(args.taxon_name)


if args.treebase:
    cites_phyl = cites_phyl + "Only returning studies with TreeBase links\n"


studies = dict()
trees = dict()
treebase_studies = set()
sys.stdout.write("Gathering references (slow)\n")
for study in phylesystem_studies_resp.response_dict['matched_studies']:
    sys.stdout.write('.')
    sys.stdout.flush()
    study_id = study['ot:studyId']
    studies[study_id] = dict()
    study_info = OT.get_study(study_id)
    nexson = study_info.response_dict['data']
    data_deposit = nexson['nexml'].get(u'^ot:dataDeposit')
    data_deposit_url = 'None'
    if data_deposit:
        data_deposit_url = data_deposit[u'@href']
    studies[study_id]['data_deposit_url'] =  data_deposit_url
    if 'treebase' in data_deposit_url:
        treebase_studies.add(study_id)
    trees[study_id] = []
    for tree in study['matched_trees']:
        treeid = tree['ot:treeId']
        trees[study_id].append(treeid)
    cites = []
    studies[study_id]['opentree_url'] = "https://tree.opentreeoflife.org/curator/study/view/{}".format(study_id)
    studies[study_id]['reference'] = nexson['nexml'].get('^ot:studyPublicationReference', 'no ref') 
    studies[study_id]['doi'] = nexson['nexml'].get('^ot:studyPublication', 'no study pub')
    if args.treebase and study_id not in treebase_studies:
            continue
    cites_phyl = cites_phyl + "\nStudy {} tree(s) {}\n".format(study_id, ', '.join(trees[study_id]))
    cites_phyl = cites_phyl + "OpenTreeUrl: " + studies[study_id]['opentree_url'] + '\n'
    cites_phyl = cites_phyl + "Reference: " + studies[study_id]['reference'] + '\n'
    cites_phyl = cites_phyl + "Data Deposit URL: " + studies[study_id]['data_deposit_url'] + '\n'



if args.output:
    ofi = open(args.output, 'w')
    ofi.write(cites_phyl)
else:
    print(cites_phyl)


