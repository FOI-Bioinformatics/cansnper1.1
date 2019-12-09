from setuptools import setup, find_packages
from os import path
import sys
import pkg_resources  # part of setuptools

#print(pkg_resources.require("CanSNPer"))
from CanSNPer import __version__
#version = __version__ #pkg_resources.require("CanSNPer")[0].version
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="CanSNPer",
    version=__version__,

    long_description=long_description,
    long_description_content_type="text/markdown",

    url="https://github.com/FOI-Bioinformatics/cansnper1.1",
     # Author details
    author='David Sundell',
    author_email='david.sundell@foi.se',

    description="CanSNPer: A toolkit for SNP-typing using NGS data.",
    license="GPL'",
    keywords="Bioinformatics SNP-typing sequence-data",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Programming Language :: Python :: 3',

        #'Intended Audience :: Developers bioinformatics',
        #'Topic :: Software Development :: Build Tools',

	"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent"
    ],
    install_requires=['numpy', 'ete3'],
    #packages=['CanSNPer'],
    #py_modules=['CanSNPer.modules.ParseXMFA',"CanSNPer.modules.DatabaseConnection"],
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    entry_points={'console_scripts': ['CanSNPer=CanSNPer.__main__:main']})
