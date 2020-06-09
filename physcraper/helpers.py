import os
import sys
import subprocess
import contextlib

if sys.version_info[0] < 3:
    str_type = unicode
else:
    str_type = bytes 


_DEBUG = 0
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)



def get_raxml_ex():
            if subprocess.check_call(["which", "raxmlHPC"]) == 0:
                rax_ex = "raxmlHPC"
            else:
                sys.stderr.write("Did not find raxml executable. Exiting \n")
                sys.exit()
            return rax_ex


def to_string(input):
    if isinstance(input, str):
        return input
    elif isinstance(input, str_type):
        output = input.decode('ascii','replace')
        return output
    else:
        return input


@contextlib.contextmanager
def cd(path):
    # print 'initially inside {0}'.format(os.getcwd())
    CWD = os.getcwd()
    os.chdir(path)
    # print 'inside {0}'.format(os.getcwd())
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        # print 'finally inside {0}'.format(os.getcwd())
        os.chdir(CWD)


#def generate_from_run(workdir,
#                      seqaln='physcraper.fas',
                      # mattype='fasta',
                      # configfi='config.out',
                      # treefile='physcraper.tre',
                      # schema_trf='newick',
                      # search_taxon = 'mrca.txt'):
                      

def write_filterblast_db(workdir, seq_name, seq, fn):
    """Writes local blast db which will be read by run_filter_blast.

    :param workdir: working directory
    :param seq_name: sequence identifier
    :param seq: sequence to write
    :param fn: file name
    :return: files with sequences written to it in fasta format
    """
    """
    """
    if not os.path.exists("{}/blast".format(workdir)):
        os.makedirs("{}/blast/".format(workdir))
    fnw = "{}/blast/{}_db".format(workdir, fn)
    fi_o = open(fnw, "a")
    fi_o.write(">{}\n".format(seq_name))
    fi_o.write("{}\n".format(str(seq).replace("-", "")))
    fi_o.close()



def standardize_label(item):
    """Make sure that the tip names are unicode.

    Function is only used if own files are used for the OtuJsonDict() function.

    :param item: original tip name
    :return: tip name in unicode
    """
    item_edit = item.replace(" ", "_")
    return item_edit



def concat(genelistdict, workdir_comb, email, num_threads=None, percentage=0.37, user_concat_fn=None, backbone=False):
    """This is to concatenate different physcraper runs into a single alignment and tree.
    genelistdict is a dict with gene names as key and the corresponding workdir
    """
    license_print()

    if not os.path.exists(path="{}/concat_checkpoint.p".format(workdir_comb)):
        if not os.path.exists(path="{}/load_single_data.p".format(workdir_comb)):
            # save_copy_code(workdir_comb)
            conc = Concat(workdir_comb, email)
            conc.concatfile = user_concat_fn
            for item in genelistdict.keys():
                conc.load_single_genes(genelistdict[item]["workdir"], genelistdict[item]["pickle"], item)
            conc.combine()
        else:
            sys.stdout.write("load single data dump file\n")
            conc = pickle.load(open("{}/load_single_data.p".format(workdir_comb), "rb"))
            # conc.dump()
        conc.sp_seq_counter()
        conc.get_largest_tre()
        conc.make_sp_gene_dict()
        conc.make_alns_dict()
        conc.concatenate_alns()
        conc.get_short_seq_from_concat(percentage)
        conc.remove_short_seq()
        conc.dump()
    else:
        sys.stdout.write("load concat_checkpoint dump file\n")
        conc = pickle.load(open("{}/concat_checkpoint.p".format(workdir_comb), "rb")) 
    conc.backbone = backbone
    conc.make_concat_table()
    conc.write_partition()
    conc.write_otu_info()
    conc.place_new_seqs(num_threads)
    
    if backbone is False:
        conc.calculate_bootstrap(num_threads)
        conc.write_labelled('RAxML_bestTree.autoMRE_fa')
    else:
        conc.est_full_tree(num_threads)
        conc.write_labelled('RAxML_bestTree.backbone_concat')
    return conc