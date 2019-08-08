import csv
import os

_DEBUG_MK = 0


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)


"""
Functions that write out additional sampling information for a PhyScraper run."""


def get_additional_GB_info(physcraper_obj):
    """Retrieves additional information given during the Genbank sequence submission
    for all included sequences and writes them out to file.

    param physcraper_obj: PhyScraper Scrape object
    """
    debug("get_additional_GB_info")
    table_keys = [
        "Genbank accession",
        "Species name",
        "authors",
        "journal",
        "publication title",
        "voucher information",
        "clone",
        "country",
        "isolate"
    ]
    with open("{}/Genbank_information_added_seq.csv".format(physcraper_obj.workdir), "w+") as output:
        writer = csv.writer(output)
        writer.writerow(table_keys)
        for entry in physcraper_obj.data.otu_dict.keys():
            # debug(entry)
            # debug(self.data.otu_dict[entry]['^physcraper:status'].split(' ')[0])
            if physcraper_obj.data.otu_dict[entry]['^physcraper:status'].split(' ')[0] not in physcraper_obj.seq_filter:
                # debug(physcraper_obj.data.otu_dict[entry].keys())
                if '^ncbi:accession' in physcraper_obj.data.otu_dict[entry].keys():
                    # debug("add info")
                    gb_id = physcraper_obj.data.otu_dict[entry]['^ncbi:accession']
                    read_handle = physcraper_obj.ids.entrez_efetch(gb_id)
                    ncbi_sp = None
                    voucher = None
                    clone = None
                    country = None
                    isolate = None
                    # debug(read_handle[0])
                    # debug(read_handle[0]["GBSeq_references"][0])
                    gb_list = read_handle[0]["GBSeq_feature-table"][0]["GBFeature_quals"]
                    # debug(gb_list)
                    for item in gb_list:
                        if item[u"GBQualifier_name"] == "organism":
                            ncbi_sp = str(item[u"GBQualifier_value"])
                            ncbi_sp = ncbi_sp.replace(" ", "_")
                        if item[u"GBQualifier_name"] == "specimen_voucher":
                            voucher = str(item[u"GBQualifier_value"])
                        if item[u"GBQualifier_name"] == "clone":
                            clone = str(item[u"GBQualifier_value"])
                        if item[u"GBQualifier_name"] == "country":
                            country = str(item[u"GBQualifier_value"])
                        if item[u"GBQualifier_name"] == "isolate":
                            isolate = str(item[u"GBQualifier_value"])
                    # debug(read_handle[0])
                    if "GBSeq_references" in read_handle[0].keys():
                        authors = read_handle[0]["GBSeq_references"][0][u'GBReference_authors']
                        journal = read_handle[0]["GBSeq_references"][0][u'GBReference_journal']
                        publication = read_handle[0]["GBSeq_references"][0][u'GBReference_title']
                        info = [gb_id, ncbi_sp, authors, journal, publication, voucher, clone, country, isolate]
                        writer.writerow(info)


def write_otu_info(physcraper_obj):
    """Writes output table to file

    1. a file with all relevant GenBank info to file (otu_dict).

    :param physcraper_obj: PhyScraper Scrape object
    :return: writes output to file
    """
    debug("write out infos")
    otu_dict_keys = [
        "^ncbi:gi",
        "^ncbi:accession",
        "^ot:originalLabel",
        "^physcraper:last_blasted",
        "^physcraper:status",
        "^physcraper:TaxonName",
        "^ncbi:title",
        "^ncbi:taxon",
        "^ncbi:TaxonName",
        "^ot:ottId",
        "^ot:ottTaxonName"
    ]
    with open("{}/otu_seq_info.csv".format(physcraper_obj.workdir), "w+") as output:
            writer = csv.writer(output)
            wr = ["otuID"]
            for key in otu_dict_keys:
                wr.append(key)
            writer.writerow(wr)
    with open("{}/otu_seq_info.csv".format(physcraper_obj.workdir), "a") as output:
        writer = csv.writer(output)
        for otu in physcraper_obj.data.otu_dict.keys():
            rowinfo = [otu]
            for item in otu_dict_keys:
                if item in physcraper_obj.data.otu_dict[otu].keys():
                    tofile = str(physcraper_obj.data.otu_dict[otu][item]).replace("_", " ")
                    rowinfo.append(tofile)
                else:
                    rowinfo.append("-")
            writer.writerow(rowinfo)


def taxon_sampling(filterblast_obj, downtorank=None):
    """Write out file which contains the taxon smapling.

    Writes output table to file: table with taxon names and sampling
    It uses the self.sp_d to get sampling information, that's why the downtorank is required.

    :param filterblast_obj: FilterBlast object
    :param downtorank: hierarchical filter
    :return: writes output to file
    """
    debug("write out taxon sampling")
    sp_d = filterblast_obj.sp_dict(downtorank)
    sp_info = {}
    for k in sp_d:
        sp_info[k] = len(sp_d[k])
    with open("{}/taxon_sampling.csv".format(filterblast_obj.workdir), "w") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in sp_info.items():
            # print(filterblast_obj.ids)
            if filterblast_obj.config.blast_loc == "remote":
                spn = filterblast_obj.ids.ncbiid_to_spn[key]
            else:
#                spn = filterblast_obj.ids.ncbiid_to_spn[key]
                spn = filterblast_obj.ids.ncbi_parser.get_name_from_id(key) #TODO this was throwing a pandas error
            writer.writerow([key, spn, value])


def write_not_added_info(physcraper_obj, gb_id, reason=None):
    """Writes out infos of not added seq based on information provided in reason.

    Is not used, as the output file can easily get to 100GB.

    :param physcraper_obj: PhyScraper Scrape object
    :param item:  retrieved seq that was not added
    :param reason: optional argument to give reason for not adding
    :return: writes output to file
    """
    debug("write not added infos")
    tab_keys = [
        "^ncbi:gi",
        "^accession",
        "sscinames",
        "staxids",
        # "title",
        "length",
        # "hsps",
        # "pident",
        "evalue",
        "bitscore"
        # "sseq"
    ]
    if not os.path.exists(path="{}/info_not_added_seq.csv".format(physcraper_obj.workdir)):
        with open("{}/info_not_added_seq.csv".format(physcraper_obj.workdir), "w+") as output:
            writer = csv.writer(output)
            writer.writerow(tab_keys)
    with open("{}/info_not_added_seq.csv".format(physcraper_obj.workdir), "a") as output:
        writer = csv.writer(output)
        if gb_id in physcraper_obj.data.gb_dict:
            rowinfo = []
            for key in tab_keys:
                tofile = str(physcraper_obj.data.gb_dict[gb_id].get(key,"-")).replace("_", " ")
                rowinfo.append(tofile)
        else:
            rowinfo = [gb_id,'-','-','-','-','-','-']
        rowinfo.append(reason)
        writer.writerow(rowinfo)


def write_not_added(ncbi_id, tax_name, gb_id, reason, workdir):
    """Writes out infos of not added seq based on information provided in reason.

    :param ncbi_id:
    :param tax_name:
    :param gb_id:
    :param reason:
    :param workdir:
    :return:
    """
    debug("write not added")
    tab_keys = [
        "ncbi_id",
        "tax_name",
        "gb_id",
        "reason"
    ]
    if not os.path.exists(path="{}/not_added_seq.csv".format(workdir)):
        with open("{}/not_added_seq.csv".format(workdir), "w+") as output:
            writer = csv.writer(output)
            writer.writerow(tab_keys)
    with open("{}/not_added_seq.csv".format(workdir), "a") as output:
        writer = csv.writer(output)
        rowinfo = [ncbi_id, tax_name, gb_id, reason]
        writer.writerow(rowinfo)
    #
    # fn = open("{}/not_added_seq.csv".format(self.workdir), "a+")
    # fn.write(
    #     "not_part_of_mrca, {}, rankid: {}, ncbi_id:{}, tax_name:{}\n".format(gb_id, input_rank_id, ncbi_id, tax_name))
    # fn.close()
