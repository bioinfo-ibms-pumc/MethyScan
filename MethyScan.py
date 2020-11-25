#########################################################################
# File Name: MethyScan.py
# > Author: CaoYinghao
#########################################################################
#! /usr/bin/python

from time import ctime
import glob
import sys
import re
import os
import argparse
from datetime import datetime
from utils import PrimerDesign as pd
from utils import Threshold
from formatdb import DBFormatter
from searchdb import Searcher

def init():
	desc = """Program: MethyScan
Version: 1.0
Contact: Yinghao Cao <yhcao@ibms.pumc.edu.cn>
	"""

	usage = "%(prog)s"
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc,usage=usage)

	subparser = parser.add_subparsers(help="command",dest="command")

	designparser = subparser.add_parser('designprimer',help="Design primer pairs for target sequence.")
	designparser.add_argument('-i','--inputfile',required=True,help="Target sequence file in fasta format.")
	designparser.add_argument('-o','--outputfile',required=True, help="Output file.")
	designparser.add_argument('--minlength',default=20,help="Primer minimum length.")
	designparser.add_argument('--maxlength',default=25,help="Primer maximum length.")
	designparser.add_argument('--mintm',default=50,help="Primer minimum TM.")
	designparser.add_argument('--maxtm',default=60,help="Primer maximum TM.")
	designparser.add_argument('-m','--maxnum',default=3,help="Maximum primer pairs for one CpG island. [Default:%(default)s]")
	designparser.add_argument('--tmdiff',default=5,help="Maximum TM difference in primer pairs.")


	dbparser = subparser.add_parser('formatdb',help="Create methylated databse.")
	dbparser.add_argument('-g','--genomefile',required=True,help="Database file in fasta format.")
	dbparser.add_argument('-o','--outputfile',help="Output database name.")
	dbparser.add_argument('-t','--threads',default=4,help="Number of threads used for building index in bowtie. [Default:%(default)s]")

	searchparser = subparser.add_parser('searchdb',help="Search primer pairs against methylated database.")
	searchparser.add_argument('-i','--inputfile',required=True,help="Input primer sequence file in tab format likes \"amp1	p5seq	p3seq\")")
	searchparser.add_argument('-d','--database',required=True,help="Methylated database filename prefix.")
	searchparser.add_argument('-o', '--outputfile', required=True, help="Output file")
	searchparser.add_argument('-c','--mismatch',default=3,help="Maximum mismatch in one primer sequence, (<=3).")
	searchparser.add_argument('--minlength',default=100,help="Minimum amplicon length considered in mapping. [default:%(default)s]")
	searchparser.add_argument('--maxlength',default=500,help="Maximum amplicon length considered in mapping. [default:%(default)s]")
	#searchparser.add_argument('-f', '--force', action="store_true",default=False, help="Reuse the old searched files. [default:%(default)s]")
	searchparser.add_argument('-k', '--searchtype',default="ALL", help="Search against methylated database with specific strand. Type: ALL,[MCT],[MGA],[UCT],[UGA]. [Default:%(default)s]")
	searchparser.add_argument('-t', '--threads', default=1,
							  help="Number of threads used for searching. [default: %(default)s]")

	return parser

def main_cmd(args):

	now = datetime.now()

	if args.command == "designprimer":
		print("{0:40}{1:>20}".format("Design primer start!", ctime()))
		#threshold = threshold(regionstart=4700,regionend=5000)
		thresh = Threshold(maxlen=args.maxlength,minlen = args.minlength,tmdiff=args.tmdiff,mintm = args.mintm,maxtm=args.maxtm)
		mypd = pd(args.inputfile,thresh)
		mypd.designprimers(args.outputfile)
	elif args.command == "formatdb":
		formatter = DBFormatter(file=args.genomefile,oname=args.outputfile,threads=args.threads)
		formatter.convertSeq()

	elif args.command == "searchdb":
		searcher = Searcher(inputfile=args.inputfile,outputfile=args.outputfile,
				dbfile=args.database,threads=args.threads,searchtype = args.searchtype,
				mismatch=args.mismatch,minlen=args.minlength,maxlen=args.maxlength)
		searcher.alignment()
		#print("{0:40}{1:>20}".format("Evaluation on database start!", ctime()))

	print("{0:40}{1:>20}".format("Job done!", ctime()))
	costs = (datetime.now() - now).seconds
	print("{0:20}{1}{2}".format("Total costs:", costs, "s"))


if __name__ == "__main__":
	parser = init()
	args = parser.parse_args()

	if args.command not in ("designprimer","formatdb","searchdb"):
		parser.print_help()
	else:
		main_cmd(args)
