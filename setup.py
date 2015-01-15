#!/usr/bin/env python

from distutils.core import setup

setup(
    name='bamcover',
    version='0.1.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['bamcover'],
    entry_points={
        'console_scripts': ['bamcover = bamcover.bamcover:main']
    },
    url='https://github.com/bjpop/bamcover',
    license='LICENSE.txt',
    description=('FIXME'),
    long_description=('FIXME'),
    install_requires=[
        "pysam >= 0.7.8",
        "bx-python >= 0.7.1"
    ],
)
