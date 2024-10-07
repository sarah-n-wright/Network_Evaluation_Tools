#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
import re
from setuptools import setup, find_packages


with open(os.path.join('neteval', '__init__.py')) as ver_file:
    for line in ver_file:
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['argparse>=1.4.0',
               'fuzzywuzzy>=0.18.0',
               'goatools>=1.3.1',
               'hidef>=1.1.5',
                'httplib2>=0.20.2',
                'matplotlib>=3.5.0',
                'mygene>=3.2.2',
                'ndex2>=3.5.0',
                'networkx>=2.6.3,<3.0',
                'numpy>=1.21.4',
                'obonet>=1.0.0',
                'pandas>=1.3.4,<2.0',
                'python>=3.10.0',
                'requests>=2.26.0',
                'scikit-learn>=1.0.1',
                'scipy>=1.7.2',
                'seaborn>=0.13.0',
                'statsmodels>=0.13.5',
                'tqdm>=4.62.3']

setup_requirements = []

test_requirements = []

setup(
    author="Sarah Wright",
    author_email='snwright@ucsd.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    description="Package for standardization and evaluation of biological networks",
    install_requires=requirements,
    python_requires='>=3.10',
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type = 'text/x-rst',
    include_package_data=True,
    keywords=['neteval', 'network', 'interaction', 'bioinformatics'],
    name='neteval',
    packages=find_packages(include=['neteval', 'neteval.*']),
    package_dir={'neteval': 'neteval'},
    scripts=[ 'neteval/netevalcmd.py', 
             'neteval/process_data.py', 
             'neteval/prepare_evaluation_data.py', 
             'neteval/run_network_evaluation.py',
             'neteval/network_constructor.py',
             'neteval/edge_prediction.py',
             'neteval/alphafold_results.py'
            ],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/sarah-n-wright/Network_Evaluation_Tools',
    version=version,
    zip_safe=False)

