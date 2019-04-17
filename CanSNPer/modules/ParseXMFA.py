#!/usr/bin/env python3 -c

'''
ParseXMFA is a script that parse xmfa files and extract all SNPs
	The script allows an option using a flanking score that limits
	SNPs near edges of a block to have a impact on classifying decisions.
'''

__version__ = "0.1.3"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2019-04-17"
__status__ = "Production"

import argparse
import os, string
from CanSNPer.modules.DatabaseConnection import DatabaseConnection

class CanSNPerClassification(object):
	"""docstring for CanSNPerClassification."""
	def __init__(self, database):
		super(CanSNPerClassification, self).__init__()
		self.database = DatabaseConnection(database)

	def snp_lister(self, organism,reference="SCHUS4.2"):
		'''Returns a list of all SNPs, their positions.

		Keyword arguments:
		organism -- the name of the organism
		return results as a dictionary with tuple SNP for each position {pos: (refBase, TargetBase)}
		'''
		res = self.database.query(
			"""SELECT Strain, Position, Derived_base, Ancestral_base, SNP
					FROM {organism}
					WHERE Strain = ?
			""".format(organism=organism), reference)
		results = {}
		for strain, pos,tbase,rbase, SNP in res.fetchall():
			results[pos] = tuple([pos,rbase, tbase,SNP])
		return results



class ParseXMFA(object):
	"""docstring for ParseXMFA."""
	def __init__(self, main=False, verbose=False, **kwargs):
		super(ParseXMFA, self).__init__()
		### Define translation table
		self.rcDict = {
			 "A" : "T",
			 "T" : "A",
			 "C" : "G",
			 "G" : "C",
			 "-" : "-",
			 "R" : "Y",
			 "Y" : "R",
			 "W" : "S",
			 "S" : "W",
			 "K" : "M",
			 "M" : "K",
			 "D" : "H",
			 "H" : "D",
			 "V" : "B",
			 "B" : "V",
			 "X" : "X",
			 "N" : "N"
		}
		'''Define all base variables'''
		self.verbose = verbose
		if main:
			self.xmfa = kwargs["xmfa"][0]
			self.database = kwargs["database"]
			self.organism = kwargs["organism"]
			self.reference = kwargs["reference"]

			'''Create connection to SNP database and retrieve registered SNPs and save first snp to look for'''
			self.SNP_DB = CanSNPerClassification(self.database)
			self.snplist, self.snp_positions = self.get_snps()
			self.current_snp = self.snp_positions.pop(0)
		## as the snps are ordered according to the reference mauve positions they can be used sequentialy without search

		'''SNPs will be stored as a sorted set containing (position, refBase, targetBase,SNPname)'''
		self.SNPS = {}

		### This function is not implemented
		'''If the score option is selected a dictionary from each set to a score is created'''
		#self.score = {}

	def get_snps(self):
		'''Retrieve snps from database object'''
		self.snplist = self.SNP_DB.snp_lister(self.organism, self.reference)
		self.snp_positions = list(self.snplist.keys())
		self.snp_positions.sort()
		return self.snplist,self.snp_positions

	def reverse_complement(self,dna):
		'''Complement and reverse DNA string'''
		dna_rev = [ self.rcDict[x] for x in dna[::-1] ]
		return "".join(dna_rev)

	def get_SNPs(self,ref,target,snp=0,head=0):
		'''Walk through the paired sequences and save SNPs'''
		snppos,rbase,tbase,_snp = snp    ## Retrieve all information known about the next upcoming SNP
		snppos -= int(head["start"])-1  ## python counts from 0 not 1
		_rbase = False                    ## create place holder for reference base
		_snp = False                    ## create place holder for target base

		'''i will count the relative position in the reference (without -)'''
		i = 0
		''' ii is the actual position in the file '''
		base_positions = range(len(ref))
		'''If sign is "-" it means the reference sequence is the reverce complement, hence positions need to be inverse'''
		if head["sign"] == "-":
			base_positions = reversed(base_positions)
		for ii in base_positions:
			if ref[ii] != "-":          ## if the reference contains a - do not count as base, else follow reference position
				i+=1
			if i == snppos:                ## check the snp position
				_snp = target[ii]        ## get base in target sequence
				if head["sign"] == "-":
					'''Again if sign is "-" the complement base needs to be retrieved'''
					_snp = self.reverse_complement(_snp)
		SNP = {snp[3]:0} ## SNP not found
		if tbase == _snp:                ## If Derived (target) base is in target sequence then call SNP
			if self.verbose: print((_rbase), (tbase), snp, head["sign"], "Called")
			SNP[snp[3]] = 1 		## Derived SNP
		elif rbase == _snp:  		## SNP is found but confirmed to be ancestral
			SNP[snp[3]] = 2			## Ancestral SNP
		return SNP

	def parse_head(self,head):
		'''This help function parses the head of an xmfa file and returns info'''
		cols = head.split(" ")
		ref = False
		sign = cols[1]
		### path = cols[2] not saved
		ref,position = cols[0].split(":")
		start, end = list(map(int,position.split("-")))
		relend = 1
		if sign == "-":
			relend *=-1
		if cols[0] == "1":
			ref = True
		return {"ref":ref,"sign":sign, "start":start,"end":end,"relend":relend}

	def read_sequence(self,seqP,res={}):
		'''read information in sequence pair'''
		seqLines = seqP.strip().split("> ")
		headinfo = seqLines[0]
		if len(seqLines) > 2:  ## Both target and reference have sequence, parse and find ----- in data
			refSeq = seqLines[1].split("\n")
			refHead = self.parse_head(refSeq.pop(0))
			targetSeq = seqLines[2].split("\n")
			targetHead = self.parse_head(targetSeq.pop(0))
			if refHead["start"] < self.current_snp and refHead["end"] > self.current_snp:
				while self.current_snp < refHead["end"]:
					ref = "".join(refSeq)
					target = "".join(targetSeq)
					offset = 1
					res = dict(**res, **self.get_SNPs(ref,target,snp=self.snplist[self.current_snp],head=refHead))
					if len(self.snp_positions) == 0:
						break
					self.current_snp = self.snp_positions.pop(0)
				return res
		return {}

	def read_xmfa(self,f=False):
		'''read xmfa file'''
		if not f:
			f = self.xmfa
		with open(f) as fin:
			fcontent = fin.read()
			#### Each aligned sequence part is separated by = sign, so start with splitting the sequences into pairs of sequence
			seqPairs = fcontent.strip().split("=")
			for seqP in seqPairs:
				### Join together all SNPs found in data
				self.SNPS = dict(**self.SNPS, **self.read_sequence(seqP))
		return self.SNPS

	def get_references(self,database):
		'''Get query'''
		db = DatabaseConnection(database)
		query = """SELECT DISTINCT(Strain) FROM Sequences"""
		res = db.query(query).fetchall()
		db.disconnect()
		return res

	def run(self, database, xmfa, organism,reference):
		'''Run script function if class is called from other script'''
		self.database = database
		self.organism = organism
		self.reference = reference
		self.xmfa = xmfa
		#print(self.database,self.organism, self.reference,self.xmfa)
		'''Create connection to SNP database and retrieve registered SNPs and save first snp to look for'''
		self.SNP_DB = CanSNPerClassification(self.database)
		self.snplist, self.snp_positions = self.get_snps()
		if len(self.snp_positions) == 0:
			print("Error no SNPs found")
			exit()
		self.current_snp = self.snp_positions.pop(0)
		return self.read_xmfa()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Parse xmfa files and extract non overlapping sequences')
	parser.add_argument('xmfa', metavar='', type=str, nargs='+', help='fasta xmfa to be parsed')
	parser.add_argument('database', metavar='',  help='CanSNP database')
	parser.add_argument('--organism', metavar='', default="Francisella", help="Specify organism")
	parser.add_argument('--reference', metavar='', default="FSC200",choices=["FSC200","SCHUS4.1","SCHUS4.2","OSU18"], help="Specify organism")
	parser.add_argument('--flank', metavar='',default=200, type=int, help="Min length of gap to be saved to output")
	parser.add_argument('--verbose',action='store_true',help="print process info, default no output")
	args = parser.parse_args()
	if args.verbose:
		print(args)
	xmfa = ParseXMFA(database=args.database, xmfa=args.xmfa, organism=args.organism,reference=args.reference, main=True)
	xmfa.read_xmfa()
