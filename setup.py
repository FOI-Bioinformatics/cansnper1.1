from setuptools import setup
import sys
import pkg_resources  # part of setuptools
version = pkg_resources.require("CanSNPer")[0].version

setup(
    name="CanSNPer",
    version=version,
    url="https://git-int.foi.se/bioinfo/cansnper1.1",
    description="CanSNPer: A toolkit for SNP-typing using NGS data.",
    license="GPL'",
    keywords="Bioinformatics SNP-typing sequence-data",
    classifiers=[
        'Development Status :: 5 - Alpha',
        'License :: OSI Approved :: GPL',
        'Programming Language :: Python :: 3'
        ],
    install_requires=['numpy', 'ete3'],
    packages=['CanSNPer'],
    py_modules=['CanSNPer.modules.ParseXMFA',"CanSNPer.modules.DatabaseConnection"],
    entry_points={'console_scripts': ['CanSNPer=CanSNPer.__main__:main']})
