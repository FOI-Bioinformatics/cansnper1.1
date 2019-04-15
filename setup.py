from setuptools import setup
import sys


setup(
    name="CanSNPer",
    version="1.1.0",
    url="https://github.com/adrlar/CanSNPer",
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
    scripts=['CanSNPer/x2fa.py'],
    entry_points={'console_scripts': ['CanSNPer=CanSNPer.__main__:main']})
