from setuptools import setup
import sys

if sys.version_info < (3,6):
    sys.exit('Sorry, Python < 3.6 is not supported')

setup(
    name='bioinfo_tools',
    version='0.2',
    url='https://github.com/sebriois/bioinfo_tools',
    author='Sebastien Briois, Guillaume Tiberi',
    author_email='sebriois@gmail.com',
    packages=[
        'bioinfo_tools',
        'bioinfo_tools.genomic_features',
        'bioinfo_tools.parsers',
        'bioinfo_tools.utils',
    ],
    keywords='bioinformatics',
    license='BSD',
    description='Python library that parses GFF, Fasta files into python classes',
    long_description=open('README.txt').read(),
    install_requires=[ "biopython" ],
    classifiers=[
        'Environment :: Web Environment',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)