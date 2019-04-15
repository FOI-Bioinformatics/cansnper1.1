#!/usr/bin/env python3 -c

'''
ParseXMFA is a script that parse xmfa files and extract all SNPs
	The script allows an option using a flanking score that limits
	SNPs near edges of a block to have a impact on classifying decisions.
'''

__version__ = "0.1.0"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se", "david.sundell@foi.se"]
__date__ = "2019-03-11"
__status__ = "Production"

'''
### Test set modifications ###
1: row 1:   base 30
2: row 5 reverse complement:   base -30  A-T
3: row 3: base 5 ->  10 bases removed
4: row 7: 10 bases added
5: row 8: SNP 3 from end T to C
'''

import argparse
import os, string
from DatabaseConnection import DatabaseConnection

import time
start_time = time.time()
current_time = start_time
def report_time(prev_time, final=False,bootstrap=1):
	current = time.time()
	diff = current - prev_time
	if bootstrap > 1:
		diff = diff/bootstrap
	seconds = diff % 60
	minutes = int((diff - seconds)/60)
	mtext = "minute"
	if minutes != 1:
		mtext +="s"
	if final:
		print("--- Time summary  {} {} {} seconds---".format(minutes,mtext, seconds))
	else:
		print("--- process finished in {} {} {} seconds---\n".format(minutes,mtext, seconds))
	return current

class CanSNPerClassification(object):
	"""docstring for CanSNPerClassification."""
	def __init__(self, database):
		super(CanSNPerClassification, self).__init__()
		self.database = DatabaseConnection(database)

	def snp_lister(self, organism,reference="SCHUS4.2"):
		'''Returns a list of all SNPs, their positions.

		Keyword arguments:
		organism -- the name of the organism

		1. FSC200
		2. SCHUS4.1
		3. SCHUS4.2
		4. OSU18
		WHERE Strain = '{reference}'
		return results as a dictionary with tuple SNP for each position {pos: (refBase, TargetBase)}
		'''
		res = self.database.query(
			"""SELECT Strain, Position, Derived_base, Ancestral_base, SNP
					FROM {organism}
					WHERE Strain = '{reference}'
			""".format(organism=organism, reference=reference))
		results = {}
		for strain, pos,tbase,rbase, SNP in res.fetchall():
			#print(strain)
			#print(tuple([pos,rbase, tbase]))
			results[pos] = tuple([pos,rbase, tbase,SNP])
		return results



class ParseXMFA(object):
	"""docstring for ParseXMFA."""
	def __init__(self, *args, **kwargs):
		super(ParseXMFA, self).__init__()
		### Define translation table
		#(pos, ref_base, t_base)
		#(pos, kmer)
		self.verbose = False
		#setA & setB
		self.xmfa = kwargs["xmfa"][0]
		self.database = kwargs["database"]
		self.organism = kwargs["organism"]
		self.SNPS = set()
		SNPs = CanSNPerClassification(self.database)
		self.snplist = SNPs.snp_lister(self.organism, kwargs["reference"])
		#print(self.snplist)
		self.snp_positions = list(self.snplist.keys())
		self.snp_positions.sort()
		self.current_snp = self.snp_positions.pop(0)
		'''SNPs will be stored as a sorted set containing (position, refBase, targetBase)'''
		#self.SNPs = set()
		'''If the score option is selected a dictionary from each set to a score is created'''
		#self.score = {}
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
		walk_file = range(len(ref))
		if head["sign"] == "-":
			# target = self.reverse_complement(target)
			# ref = self.reverse_complement(ref)
			walk_file = reversed(walk_file)
		for ii in walk_file:
			if ref[ii] != "-":          ## if the reference contains a - do not count as base, else follow reference position
				i+=1
			if i == snppos:                ## check the snp position
				_rbase = ref[ii]        ## get base in reference sequence
				_snp = target[ii]        ## get base in target sequence
				if head["sign"] == "-":
					_snp = self.reverse_complement(_snp)
		if self.verbose:
			if _rbase != _snp:
				print("SNP")
		if tbase in _snp:                ## If Derived (target) base is in target sequence then call SNP
			if self.verbose: print((_rbase), (tbase), snp, head["sign"], "Called")
			return set(snp)
		elif rbase in _snp:                ## If Ancestral base in reference then Ancestral base was found, do not call
			#print(_snp, snp)
			if self.verbose: print((_rbase), (rbase),snp, head["sign"], "Ancestral base found")
		else:                            ## Other base in position
			if self.verbose: print("Other base found", snp, _snp)
		#exit()
		return set()

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

	def read_sequence(self,seqP,res=set()):
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
					res |= self.get_SNPs(ref,target,snp=self.snplist[self.current_snp],head=refHead)
					if len(self.snp_positions) == 0:
						break

					self.current_snp = self.snp_positions.pop(0)
				return res
		return set()

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
				self.SNPS |= self.read_sequence(seqP)
		return

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
	xmfa = ParseXMFA(database=args.database, xmfa=args.xmfa, organism=args.organism,reference=args.reference)
	if not args.verbose:
		xmfa.read_xmfa()
	else:
		boot = 100
		for x in range(boot):
			xmfa.read_xmfa()
		current_time = report_time(current_time, bootstrap=boot)
