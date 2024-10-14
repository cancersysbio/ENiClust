import os
import sys
#import importlib
#import importlib_resources

"""
Provides access to example files
"""

def get_data_dir():
    """
    Returns the data directory that contains reference files.
    """
    return os.path.join(os.path.dirname(__file__), 'data')

    #with importlib.resources.path('eniclust', 'data') as data_path:
    #return importlib_resources.files('eniclust').joinpath('data')


def example_filename(fn,sub_dir=None):
    """
    Return a reference file from data directory.
    This code is adapted from https://github.com/daler/pybedtools
    """
    #print(data_dir())
    #sys.exit()
    if sub_dir:
      fn = os.path.join(data_dir(), sub_dir, fn)
    else:
      fn = os.path.join(data_dir(), fn)
    #print(fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


def create_dir(dir_path):
    '''
    Create a output directory if it's not exists.
    '''

    if not os.path.exists(dir_path):
        try:
            os.makedirs(dir_path)
        except:
            sys.exit( "Output directory (%s) could not be created." % dir_path )
    return dir_path

