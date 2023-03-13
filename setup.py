#!/usr/bin/env python
from distutils.core import setup
import os
import sys

import pip
try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

def get_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "VERSION")
    version = "N/A"

    if os.path.exists(version_file):
        try:
            with open(version_file) as fp:
                version = fp.readline().strip()
        except IOError:
            pass
    return version

# At the time of writing there is an open issue on pip > 6.0
#    Where session is required parameter. Breaks backwards compatibility.
if int(pip.__version__.split(".")[0]) >= 6:
    install_reqs = parse_requirements('requirements.txt', session=False)
else:
    install_reqs = parse_requirements('requirements.txt')


install_requires = [str(ir.requirement) for ir in install_reqs]

setup(name='snapperdb',
      version=get_version(),
      description='SNP calling pipeline tools.',
      author='Tim Dallman',
      author_email='t.j.dallman@uu.nl',
      url='https://github.com/tjdallman/snapperdb',
      download_url='https://github.com/phe-bioinformatics/snapperdb/archive/v1.0.6.tar.gz',
      packages=['snapperdb','snapperdb.snpdb','snapperdb.gbru_vcf'],
      include_package_data=True,
      package_data={'snpdb': ['snapperdb/snpdb/template_snapperdb_denovo_refs_sql']},
      scripts=['snapperdb.py'],
      install_requires=install_requires
     )