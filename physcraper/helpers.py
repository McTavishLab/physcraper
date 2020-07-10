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
                      


def standardize_label(item):
    """Make sure that the tip names are unicode.

    Function is only used if own files are used for the OtuJsonDict() function.

    :param item: original tip name
    :return: tip name in unicode
    """
    item_edit = item.replace(" ", "_")
    return item_edit


