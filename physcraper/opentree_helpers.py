"""Linker Functions to get data from OpenTree"""

import json
import copy
import sys
import os
import xml
from urllib.error import HTTPError

import requests

from dendropy import DataSet, datamodel
from opentree import OT, object_conversion
from nexson.syntax import get_subtree_otus


import physcraper
from physcraper.helpers import standardize_label, to_string



_VERBOSE = 0

def set_verbose():
    """Set output verbosity
    """
    global _VERBOSE # pylint: disable=global-statement
    _VERBOSE = 1

_DEBUG = 1

def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)

physcraper_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))





phylesystemref = "McTavish EJ, Hinchliff CE, Allman JF, Brown JW, Cranston KA, Holder MT,\
                  Phylesystem: a gitbased data store for community curated phylogenetic estimates. \
                  Bioinformatics. 2015 31 2794-800. doi: 10.1093/bioinformatics/btv276\n"
synthref = "Redelings BD, Holder MT. A supertree pipeline for summarizing phylogenetic and \
            taxonomic information for millions of species. PeerJ. 2017;5:e3058. https://doi.org/10.7717/peerj.3058 \n"

def root_tree_from_synth(tree, otu_dict, base='ott'):
    # Disable all the too many locals and branches for this function
    # pylint: disable=too-many-locals, too-many-branches
    """Uses information from OpenTree of Life to suggest root.
    :param tree: dendropy :class:`Tree
    :param otu_dict: a dictionary of tip label metadata, inculding an '^ot:ottId'attribute
    'param base: either `synth` or `ott.` If `synth` will use OpenTree synthetic tree relationships to root input tree,
    if `ott` will use OpenTree taxonomy.
    """
    leaves = [leaf.taxon.label for leaf in tree.leaf_nodes()]
    spp = {otu_dict[otu]['^ot:ottId'] for otu in leaves}
    if None in spp:
        spp.remove(None)
    assert(base in ['synth', 'ott'])
    if base == 'synth':
    ## ONLY included SPP with phylo information.
        synth_ids = ottids_in_synth()
        synth_spp = set()
        for sp in spp:
            if sp in synth_ids:
                synth_spp.add(int(sp))
        if len(synth_spp) <= 3:
            sys.stdout.write("Didn't find enough taxon matches in synth tree to root. Tree is unrooted\n")
            return tree
        resp = OT.synth_induced_tree(ott_ids=synth_spp)
        induced_tree_of_taxa = resp.tree
        base = "OpenTree synthetic tree"
    elif base == 'ott':
        tax_mrca = OT.taxon_mrca(spp).response_dict['mrca']['ott_id']
        resp = OT.taxon_subtree(tax_mrca)
        taxleaves = [leaf.taxon.label for leaf in resp.tree.leaf_nodes()]
        matches = [label for label in taxleaves if int(label.split()[-1].strip('ott')) in spp]
        induced_tree_of_taxa = resp.tree.extract_tree_with_taxa_labels(matches)
        base = "OpenTree taxonomy"
    node = None
    for node in induced_tree_of_taxa:
        if node.parent_node is None:
            break # stop when node has no parents, i.e., it is the root node
    children = node.child_nodes()
    representative_taxa = []
    for child in children:
        if len(child.child_nodes()) > 1:
            subtre = child.extract_subtree()
            tip = subtre.leaf_nodes()[0]
            representative_taxa.append(tip.taxon.label)
        else:
            representative_taxa.append(child.leaf_nodes()[0].taxon.label)

    sys.stdout.write("Rooting tree based on taxon relationships in {}. \
                      \nRoot will be MRCA of {}\n".format(base, ", ".join(representative_taxa)))

    ## Get tips for those taxa:
    ## TODO: make following code a function of its own
    # def get_tips(tree, otu_dict, leaves):
    phyloref = set()
    for tax in representative_taxa:
        ott_id = int(tax.split()[-1].strip('ott_id'))
        phyloref.add(ott_id)

    tips_for_root = set()
    for otu in otu_dict:
        if otu_dict[otu]['^ot:ottId'] in phyloref:
            if  otu in leaves:
                tips_for_root.add(otu)
                continue
    mrca = tree.mrca(taxon_labels=tips_for_root)
    tree.reroot_at_node(mrca)
    return tree

def ottids_in_synth(synthfile=None):
    """Checks if ottids are present in current synth tree,
    using a file listing all current otts in synth (v12.3)
    :param synthfile: defaults to taxonomy/ottids_in_synth.txt
    """
    if synthfile is None:
        synthfile = open("{}/taxonomy/ottids_in_synth.txt".format(physcraper_dir))
    synth_ottids = set()
    for lin in synthfile:
        ottid = lin.lstrip('ott').strip()
        if len(ottid) >= 1:
            synth_ottids.add(int(ottid))
    return synth_ottids


def check_if_ottid_in_synth(ottid):
    """Web call to check if ott id in synth tree. NOT USED.
    """
    url = 'https://api.opentreeoflife.org/v3/tree_of_life/node_info'
    payload = json.dumps({"ott_id":int(ottid)})
    headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
    try:
        r = requests.post(url, data=payload, headers=headers)
        if r.status_code == 200:
            return 1
        elif r.status_code == 400:
            return 0
        elif r.status_code == 502:
            sys.stderr.write("Bad OpenTree taxon ID: {}".format(ottid))
            return 0
        sys.stderr.write("unexpected status code from node_info call: {}".format(r.status_code))
        return 0
    except requests.ConnectionError:
        sys.stderr.write("Connection Error - coud not get taxon information from OpenTree\n")

def get_ottid_from_gbifid(gbif_id):
    """Returns a dictionary mapping gbif_ids to ott_ids.
    ott_id is set to 'None' if the gbif id is not found in the Open Tree Txanomy
    """
    url = 'https://api.opentreeoflife.org/v3/taxonomy/taxon_info'
    headers = {'content-type':'application/json'}
    tax = int(gbif_id)
    payload = json.dumps(dict(source_id='gbif:{}'.format(tax)))
    res = requests.post(url, data=payload, headers=headers)
    if res.status_code == 200:
        ott_id = int(res.json()['ott_id'])
        return ott_id
    sys.stderr.write("error getting ott_id for gbif id {}, {}, {}".format(tax, res.status_code, res.reason))
    return None



def bulk_tnrs_load(filename):
    """Read in outputs from OpenTree Bulk TNRS,
    translates to a Physcraper otu_dictionary.
    :param filename: input json file
    """
    otu_dict = {}
    with open(filename) as data_file:
        input_dict = json.load(data_file)
    if "names" in input_dict.keys():
        for name in input_dict["names"]:
            i = 1
            otu = "otu" + name['id'].strip('name')
            while otu in otu_dict.keys():
                otu = "{}_{}".format(otu, i)
                i += 1
            otu_dict[otu] = {"^ot:originalLabel":name["originalLabel"]}
            if name.get("ottTaxonName"):
                otu_dict[otu]["^ot:ottTaxonName"] = name["ottTaxonName"]
            if name.get("ottId"):
                otu_dict[otu]["^ot:ottId"] = name["ottId"]
            for source in name.get("taxonomicSources", []):
                #debug(source)
                if source:
                    taxsrc = source.split(":")
                    assert len(taxsrc) == 2, taxsrc
                    otu_dict[otu]["^{}:taxon".format(taxsrc[0])] = taxsrc[1]
        for otu in otu_dict:
            otu_dict[otu]["^physcraper:status"] = "original"
            otu_dict[otu]["^physcraper:last_blasted"] = None
            otu_dict[otu]["^physcraper:ingroup"] = "unknown"
    else:
        for otu in input_dict:
            assert input_dict[otu]["^physcraper:status"], otu_dict[otu]
            otu_dict = input_dict
    return otu_dict



#def get_cite_for_study(study_id):
# to complement the function we should do also def get_ott_id_data(ott_id):

# following function assumes that study_id is an object from
# url = 'https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree'
# headers = {'content-type':'application/json'}
# payload = json.dumps(dict(ott_ids=ott_ids, label_format = label_format))
# res = requests.post(url, data=payload, headers=headers)

def get_citations_from_json(synth_response, citations_file):
    """Get ciattions for studies in an induced synth tree repsonse.
    :param synth_response: Webservice call record
    :param citations_file: output file
    """
    assert isinstance(citations_file, str)
    f = open(citations_file, "w+")
    sys.stdout.write("Gathering citations ...")
    assert 'supporting_studies' in synth_response.keys(), synth_response.keys()
    for study in synth_response['supporting_studies']:
        study = study.split('@')[0]
        index_url = 'https://api.opentreeoflife.org/v3/studies/find_studies'
        headers = {'content-type':'application/json'}
        payload = json.dumps({"property":"ot:studyId", "value":study, "verbose":"true"})
        res_cites = requests.post(index_url, data=payload, headers=headers)
        new_cite = res_cites.json()['matched_studies']
#        debug(new_cite)
        sys.stdout.write('.')
        if new_cite:
            f.write(to_string(new_cite[0].get('ot:studyPublicationReference', '')) + '\n'
                    + new_cite[0].get('ot:studyPublication', '') + '\n')
    f.close()
    sys.stdout.write("Citations printed to {}\n".format(citations_file))

# another way to do it is calling each id
# get_citation_for_study(study_id)
# use append

def conflict_tree(inputtree, otu_dict):
    """Write out a tree with labels that work for the OPenTree Conflict API
    """
    tmp_tree = copy.deepcopy(inputtree)
    i = 1
    for node in tmp_tree:
        i += 1
        if node.taxon:
            otu = otu_dict[node.taxon.label]
            ottid = otu['^ot:ottId']
            new_label = "_nd{}_ott{}".format(i, ottid)
            node.taxon.label = new_label
        else:
            node.label = "_nd{}_".format(i)
    return tmp_tree

def get_tree_from_synth(ott_ids, label_format="name", citation="cites.txt"):
    """Wrapper for OT.synth_induced_tree that also pulls citations"""
    synth_json = OT.synth_induced_tree(ott_ids=ott_ids, label_format=label_format)
    get_citations_from_json(synth_json.response_dict, citation)
    return synth_json.tree

def get_tree_from_study(study_id, tree_id, label_format="ot:originallabel"):
    """Create a dendropy Tree object from OpenTree data.
    :param study_id: OpenTree Study Id
    :param tree_id: OpenTree tree id
    :param label_format: One of 'id', 'name', "ot:originallabel", "ot:ottid", "ot:otttaxonname".
    defaults to "ot:originallabel"
    """
    assert label_format in ['id', 'name', "ot:originallabel", "ot:ottid", "ot:otttaxonname"]
    study = OT.get_study(study_id)
    study_nexson = study.response_dict['data']
    DC = object_conversion.DendropyConvert()
    tree_obj = DC.tree_from_nexson(study_nexson, tree_id, label_format)
    cites = study_nexson['nexml']['^ot:studyPublicationReference']
    return tree_obj, cites

# ATT is a dumb acronym for Alignment Tree Taxa object
def generate_ATT_from_phylesystem(alnfile,
                                  aln_schema,
                                  workdir,
                                  configfile,
                                  study_id,
                                  tree_id,
                                  search_taxon=None,
                                  tip_label='^ot:originalLabel'):
    # Disable all the too many locals, branches and statements for this function
    # pylint: disable=too-many-locals, too-many-branches, too-many-statements
    """
    Gathers together tree, alignment, and study info; forces names to OTT ids.

    Study and tree ID's can be obtained by using python ./scripts/find_trees.py LINEAGE_NAME

    Spaces vs underscores kept being an issue, so all spaces are coerced to underscores when data are read in.

    :param aln: dendropy :class:`DnaCharacterMatrix
    <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix>` alignment object
    :param workdir: path to working directory
    :param config_obj: config class containing the settings
    :param study_id: OpenTree study id of the phylogeny to update
    :param tree_id: OpenTree tree id of the phylogeny to update, some studies have several phylogenies
    :param phylesystem_loc: access the GitHub version of the OpenTree
    data store, or a local clone
    :param search_taxon: optional.  OTT id of the MRCA of the clade that shall be updated
    :return: object of class ATT
    """
    assert(tip_label in ['^ot:originalLabel', 'otu', "^ot:ottTaxonName", "^ot:ottId"])
    try:
        study = OT.get_study(study_id)
    except AttributeError:
        sys.stderr.write("Failure getting study {} from OpenTree phylesystem.".format(study_id))
        sys.exit()
    try:
        study_nexson = study.response_dict['data']
        DC = object_conversion.DendropyConvert()
        tree_obj = DC.tree_from_nexson(study_nexson, tree_id)
    except KeyError:
        sys.stderr.write("Failure getting tree {} from study {} from OpenTree phylesystem.".format(tree_id, study_id))
        sys.exit()
    # this gets the taxa that are in the subtree with all of their info - ott_id, original name,
    otu_dict = {tn.taxon.otu:{} for tn in tree_obj.leaf_node_iter()}
    orig_lab_to_otu = {}
    treed_taxa = {}
    ingroup_otus = get_subtree_otus(study_nexson,
                                    tree_id=tree_id,
                                    subtree_id="ingroup",
                                    return_format="otu_id")
    if ingroup_otus is None:
        sys.stdout.write("No ingroup annotation found in tree; using all taxa.\n \
                          Please update tree annotation through OpenTree curation app.\n")
    for leaf in tree_obj.leaf_node_iter():
        tn = leaf.taxon
        otu_id = tn.otu
        otu_dict[otu_id]["^ot:ottId"] = tn.ott_id
        otu_dict[otu_id]["^ot:ottTaxonName"] = tn.ott_taxon_name
        otu_dict[otu_id]["^ot:originalLabel"] = tn.original_label.replace(" ", "_")
        otu_dict[otu_id]["^physcraper:status"] = "original"
        otu_dict[otu_id]["^physcraper:last_blasted"] = None
        if ingroup_otus:
            if otu_id in set(ingroup_otus):
                otu_dict[otu_id]["^physcraper:ingroup"] = True
            else:
                otu_dict[otu_id]["^physcraper:ingroup"] = False
        else:
            otu_dict[otu_id]["^physcraper:ingroup"] = 'unknown'
        orig = otu_dict[otu_id].get(u"^ot:originalLabel").replace(" ", "_")
        orig_lab_to_otu[orig] = otu_id
        if tip_label == 'otu':
            tn.label = otu_id
        else:
            tn.label = otu_dict[otu_id].get(tip_label)
        treed_taxa[orig] = otu_dict[otu_id].get(u"^ot:ottId")
    # need to prune tree to seqs and seqs to tree...
    ott_mrca = None
    if search_taxon:
        if search_taxon.isinstance(list):
            ott_ids = set(search_taxon)
            ott_mrca = get_mrca_ott(ott_ids)
        else:
            ott_mrca = int(search_taxon)
    if ott_mrca is None:
        ingroup_ott_ids = set()
        for otu_id in otu_dict:
            if ingroup_otus:
                if otu_dict[otu_id]["^physcraper:ingroup"]:
                    ingroup_ott_ids.add(otu_dict[otu_id].get(u"^ot:ottId"))
            else:
                ingroup_ott_ids.add(otu_dict[otu_id].get(u"^ot:ottId"))
        if None in ingroup_ott_ids:
            ingroup_ott_ids.remove(None)
        assert len(ingroup_ott_ids) >= 1
        ott_mrca = get_mrca_ott(ingroup_ott_ids)
    otu_newick = tree_obj.as_string(schema="newick")
    return physcraper.aligntreetax.AlignTreeTax(tree=otu_newick,
                                                otu_dict=otu_dict,
                                                alignment=alnfile,
                                                aln_schema=aln_schema,
                                                search_taxon=ott_mrca,
                                                workdir=workdir,
                                                configfile=configfile)
    # newick should be bare, but alignment should be DNACharacterMatrix



def get_dataset_from_treebase(study_id):
    """
        Given a tree in OpenTree with mapped tips, this function gets the corresponding
        alignment from treeBASE if available.
    """
    try:
        study = OT.get_study(study_id)
        nexson = study.response_dict['data']
    except HTTPError as err:
        sys.stderr.write(err)
        sys.stderr.write("Could not find study id {} in OpenTree phylesystem.\n".format(study_id))
    treebase_url = nexson['nexml'][u'^ot:dataDeposit'][u'@href']
    if 'treebase' not in nexson['nexml'][u'^ot:dataDeposit'][u'@href']:
        sys.stderr.write("No treeBASE record associated with study ")
        sys.exit(-2)
    else:
        tb_id = treebase_url.split(':S')[1]
        try:
            url = "https://raw.githubusercontent.com/TreeBASE/supertreebase/master/data/treebase/S{}.xml".format(tb_id)
            try:
                dna = DataSet.get(url=url, schema="nexml")
            except xml.etree.ElementTree.ParseError:
                sys.stdout.write("Error reading nexml, from supertreebase, will check TreeBASE\n")
                url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexml".format(tb_id)
                dna = DataSet.get(url=url, schema="nexml")
                return dna
        except HTTPError as err:
            try:
                url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexml".format(tb_id)
                dna = DataSet.get(url=url, schema="nexml")
            except HTTPError:
                sys.stderr.write("Alignment not found on treeBASE or supertreebase.\n")
                sys.exit()
        if _DEBUG:
            sys.stderr.write(url + "\n")
        return dna

def count_match_tree_to_aln(tree, dataset):
    """Assess how many taxa mantch between multiple genes in an alignment data set and input tree"""
    aln_match = {}
    i = 0
    leaves = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    for mat in dataset.char_matrices:
        aln_match[i] = 0
        if isinstance(mat, datamodel.charmatrixmodel.DnaCharacterMatrix):
            for tax in mat:
                if tax.label in leaves:
                    aln_match[i] += 1
        i += 1
    return aln_match

def get_max_match_aln(tree, dataset, min_match=3):
    """Select an alignment from a DNA dataset
    """
    aln_match = count_match_tree_to_aln(tree, dataset)
    max_val = min_match
    max_match = None
    for aln in aln_match:
        if aln_match[aln] > max_val:
            max_match = aln
            max_val = aln_match[aln]
    if max_match is not None:
        assert isinstance(dataset.char_matrices[max_match], datamodel.charmatrixmodel.DnaCharacterMatrix)
        return dataset.char_matrices[max_match]
    return None

def deconcatenate_aln(aln_obj, filename, direc):
    """Split out seperate concatended alignments. NOT TESTED
    """
    #dna1 = dendropy.DnaCharacterMatrix.get(file=open("treebase_alns/M4358.nex"), schema="nexus")
    for label in aln_obj.character_subsets.keys():
        sys.stdout.write("deconcatenating {}".format(label))
        submat = aln_obj.export_character_subset(label)
        submat.write(path="{}/{}_{}.fasta".format(direc, filename, label), schema="fasta")


def scraper_from_opentree(study_id, tree_id, alnfile, workdir, aln_schema, configfile=None):
    """Pull tree from OpenTree to create a physcraper object
    """
    # Read in the configuration information
    data_obj = generate_ATT_from_phylesystem(alnfile=alnfile,
                                             aln_schema=aln_schema,
                                             workdir=workdir,
                                             configfile=configfile,
                                             study_id=study_id,
                                             tree_id=tree_id)
    ids = physcraper.IdDicts(data_obj.config)
    scraper = physcraper.PhyscraperScrape(data_obj, ids)
    return scraper



def OtuJsonDict(id_to_spn, id_dict):
    # Disable all the too many locals for this function
    # pylint: disable=too-many-locals
    """Makes otu json dict, which is also produced within the openTreeLife-query.

     This function is used, if files that shall be updated are not part of the OpenTreeofLife project.
    It reads in the file that contains the tip names and the corresponding species names.
    It then tries to get the different identifier from the OToL project or if not from ncbi.

    Reads input file into the var sp_info_dict, translates using an IdDict object
    using web to call Open tree, then ncbi if not found.

    :param id_to_spn: user file, that contains tip name and corresponding sp name for input files.
    :param id_dict: uses the id_dict generated earlier
    :return: dictionary with key: "otu_tiplabel" and value is another dict with the keys '^ncbi:taxon',
                                                    '^ot:ottTaxonName', '^ot:ottId', '^ot:originalLabel',
                                                    '^user:TaxonName', '^physcraper:status', '^physcraper:last_blasted'
    """
    sys.stdout.write("Set up OtuJsonDict \n")
    sp_info_dict = {}
    nosp = []
    with open(id_to_spn, mode="r") as infile:
        for lin in infile:
            ottid, ottname, ncbiid = None, None, None
            tipname, species = lin.strip().split(",")
            clean_lab = standardize_label(tipname)
            assert clean_lab not in sp_info_dict, ("standardized label ('{}') of \
                                                    `{}` already exists".format(clean_lab, tipname))
            otu_id = clean_lab
            spn = species.replace("_", " ")
            info = get_ott_taxon_info(spn)
            if info:
                ottid, ottname, ncbiid = info
            if not info:
                sys.stderr.write("match to taxon {} not found in open tree taxonomy or NCBI. "
                                 "Proceeding without taxon info\n".format(spn))
                nosp.append(spn)
            ncbi_spn = None
            if ncbiid is not None:
                ncbi_spn = spn
            else:
                ncbi_spn = id_dict.ott_to_ncbi[ottid]
            sp_info_dict[otu_id] = {
                "^ncbi:taxon": int(ncbiid),
                "^ot:ottTaxonName": ottname,
                "^ot:ottId": ottid,
                "^ot:originalLabel": tipname,
                "^user:TaxonName": species,
                "^physcraper:status": "original",
                "^physcraper:last_blasted": None,
                }
            if ncbi_spn is not None:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = ncbi_spn
                sp_info_dict[otu_id]["^ncbi:TaxonName"] = ncbi_spn
            elif ottname is not None:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = ottname
            elif sp_info_dict[otu_id]['^user:TaxonName']:
                sp_info_dict[otu_id]["^physcraper:TaxonName"] = sp_info_dict[otu_id]['^user:TaxonName']
            assert sp_info_dict[otu_id]["^physcraper:TaxonName"]  # is not None
    return sp_info_dict


#####################################
def get_nexson(study_id):
    """Grabs nexson from phylesystem"""
    study = OT.get_study(study_id)
    nexson = study.response_dict['data']
    return nexson


def get_mrca_ott(ott_ids):
    """finds the mrca of the taxa in the ingroup of the original
    tree. The blast search later is limited to descendants of this
    mrca according to the ncbi taxonomy

    Only used in the functions that generate the ATT object.

    :param ott_ids: list of all OToL identifiers for tip labels in phylogeny
    :return: OToL identifier of most recent common ancestor or ott_ids
    """
    debug("get_mrca_ott")
    mrca_node = OT.synth_mrca(ott_ids=ott_ids).response_dict
    if u'nearest_taxon' in mrca_node.keys():
        tax_id = mrca_node[u'nearest_taxon'].get(u'ott_id')
        if _VERBOSE:
            sys.stdout.write('(v3) MRCA of sampled taxa is {}\n'.format(mrca_node[u'nearest_taxon'][u'name']))
    elif u'taxon' in mrca_node['mrca'].keys():
        tax_id = mrca_node['mrca'][u'taxon'][u'ott_id']
        if _VERBOSE:
            sys.stdout.write('(v3) MRCA of sampled taxa is {}\n'.format(mrca_node['mrca'][u'taxon'][u'name']))
    else:
        sys.stderr.write('(v3) MRCA of sampled taxa not found. Please find and input an '
                         'appropriate OTT id as ingroup mrca in generate_ATT_from_files')
        sys.exit(-4)
    return tax_id

def get_ott_taxon_info(spp_name):
    """get ottid, taxon name, and ncbid (if present) from Open Tree Taxonomy.
    ONLY works with version 3 of Open tree APIs

    :param spp_name: species name
    :return:
    """
    #This is only used to write out the opentree info file. Could use NCBI id's instead of name, and likely be quicker.
    # debug(spp_name)
    try:
        call = OT.tnrs_match([spp_name], do_approximate_matching=True)
        res = call.response_dict['results'][0]
    except IndexError:
        sys.stderr.write("match to taxon {} not found in open tree taxonomy\n".format(spp_name))
        return None, None, None
    if res['matches'][0]['is_approximate_match'] == 1:
        sys.stderr.write("""exact match to taxon {} not found in open tree taxonomy.
                          Check spelling. Maybe {}?\n""".format(spp_name, res['matches'][0][u'ot:ottTaxonName']))
        return None, None, None
    if res["matches"][0]["is_approximate_match"] == 0:
        ottid = res["matches"][0]["taxon"][u"ott_id"]
        ottname = res["matches"][0]["taxon"][u"unique_name"]
        ncbi_id = None
        for source in res["matches"][0]["taxon"][u"tax_sources"]:
            if source.startswith("ncbi"):
                ncbi_id = source.split(":")[1]
        return ottid, ottname, ncbi_id
    sys.stderr.write("match to taxon {} not found in open tree taxonomy".format(spp_name))
    return None, None, None
