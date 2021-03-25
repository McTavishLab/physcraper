"""Some minor handy functions"""
import os
import sys
import subprocess
import contextlib



_DEBUG = 0
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)


def get_raxml_ex():
    """Check location of RaxML exectable"""
    if subprocess.check_call(["which", "raxmlHPC"]) == 0:
        rax_ex = "raxmlHPC"
        return rax_ex
    sys.stderr.write("Did not find raxml executable. Exiting \n")
    sys.exit()



def to_string(inputstr):
    """Coerce to string"""
    if isinstance(inputstr, bytes):
        output = inputstr.decode('ascii', 'replace')
        return output
    return inputstr


@contextlib.contextmanager
def cd(path):
    """
        Change directories and return to original directory.
        Code taken from
        https://kitchingroup.cheme.cmu.edu/blog/2013/06/16/Automatic-temporary-directory-changing/
        Mainly used for RAxML runs.
        If path does not exist, it will silently stay in the same path.
    """
    # print 'initially inside {0}'.format(os.getcwd())
    curr = os.getcwd()
    os.chdir(path)
    # print 'inside {0}'.format(os.getcwd())
    try:
        yield
    except TypeError:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        # print('Changed dir to {}'.format(os.getcwd()))
        os.chdir(curr)

def standardize_label(item):
    """Make sure that the tip names are unicode.

    Function is only used if own files are used for the OtuJsonDict() function.

    :param item: original tip name
    :return: tip name in unicode
    """
    item_edit = item.replace(" ", "_")
    return item_edit
