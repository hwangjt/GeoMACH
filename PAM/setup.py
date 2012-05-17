import os.path
import setuptools
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

sources = [
    'PAM/components/component.py', 
    'PAM/components/halfbody.py', 
    'PAM/components/fullplate.py', 
    'PAM/components/fullplate.py', 
    'PAM/configurations/configuration.py', 
    'PAM/configurations/fuse.py', 
    'PAM/configurations/wing.py', 
    'PAM/configurations/wingbody.py', 
    'PAM/configurations/wingbodytail.py',
    ]

config = Configuration(name='PAM')

kwds = {'install_requires':['numpy','scipy', 'PUBS'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'public domain',
        }
kwds.update(config.todict())

setup(**kwds)
