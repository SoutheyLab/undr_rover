#!/usr/bin/env python

from distutils.core import setup

setup(
    name='undr_rover',
    version='1.1.0',
    author='Roger Li',
    author_email='r.li3@student.unimelb.edu.au',
    packages=['undr_rover'],
    entry_points={
        'console_scripts':['undr_rover = undr_rover.undr_rover:main']
    },
    url='https://github.com/rli99/undr_rover',
    license='LICENSE.txt',
    description=(
        'UNDR ROVER: Unmapped primer directed read overlap variant caller.'),
    long_description=open('README.md').read(),
    install_requires=[
        "PyVCF==0.6.7",
        "biopython==1.66",
        "pyfaidx==0.4.7.1",
        "pysam==0.9.0"
    ],
)
