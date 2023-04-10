#!/usr/bin/env python

import sys
if sys.version_info > (3, ) and sys.version_info < (3, 7):
    sys.exit("ERROR: LRphase requires Python 3.7 or greater")

#from __future__ import absolute_import
#from __future__ import print_function

from glob import glob
from os.path import basename, splitext
from pathlib import Path

try:
    from setuptools import setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

    
def readme():
    if sys.version_info > (3, ):
        with open(Path(__file__).parent.resolve() / 'README.md', encoding='utf-8') as md:
            return md.read()
    else:
        with open('README.md') as md:
            return md.read()


def main():

    metadata = dict(
        name = 'LRphase',
        version = '1.1.1',
        license = 'MIT',
        description = 'Phasing individual long reads using known haplotype information.',
        description_content_type = 'text/plain',
        long_description = readme(),
        long_description_content_type = 'text/markdown',
        author = 'Greg Farnum',
        author_email = 'gregfar@umich.edu',
        url = 'https://github.com/Boyle-Lab/LRphase.git',
        packages = find_packages('src'),
        package_dir = {'':'src'},
        py_modules = [splitext(basename(path))[0] for path in glob('src/LRphase/*.py')],
        include_package_data = True,
        zip_safe = False,
        classifiers = [
            # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: Implementation :: CPython',
            'Programming Language :: Python :: Implementation :: PyPy',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: Microsoft',
            'Operating System :: MacOS',
            'Operating System :: Unix',
            'Programming Language :: Python',
            'Topic :: Documentation',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Software Development :: Libraries :: Python Modules',
            'Topic :: Software Development :: Version Control :: Git'
        ],
        project_urls = {
            'Issue Tracker':'https://github.com/Boyle-Lab/LRphase/issues',
        },
        keywords = [
            'long-reads', 'phasing', 'haplotype',
        ],
        python_requires = '>=3.7',
        install_requires = [
            'pysam>=0.16.0.1',
            'biopython>=1.78',
            'pyliftover>=0.4',
            'powerlaw>=1.4.6',
            'numpy>=1.20.1'
            'requests>=2.26.0'
        ],
        extras_require = {
            'mappy':['mappy'],
        },
        setup_requires=[
            'pytest-runner',        
        ],
        entry_points = {
            'console_scripts':[
            'LRphase = LRphase.cli:main',
            ]
        }
    )

    setup(**metadata)


if __name__ == "__main__":
    main()
