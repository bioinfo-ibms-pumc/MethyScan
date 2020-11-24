#########################################################################
# File Name: formatdb.py
# > Author: CaoYinghao
# > Mail: caoyinghao@gmail.com 
# Created Time: Sun Nov  1 12:18:17 2020
#########################################################################
#! /usr/bin/python

import glob
import sys
import re
import os
from datetime import datetime
from time import ctime
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess 

class DBFormatter(object):
	TYPES = ['MCT','MGA','UCT','UGA']

	def __init__(self,file,dir="./",oname = "",threads=4):

		self.gfile = file
		self.dir = dir
		self.slices = -1
		self.reflists = {}
		if oname is None or oname == "":
			oname = file.split("/")[-1]
		self.oname = oname
		self.threads = threads



	def convertSeq(self):
		#self.translate()
		self.indexdb()

	def translate(self):
		if not os.path.exists(self.gfile):
			print("Error: " + self.gfile + " is not exist!")
			exit()
		print("{0:40}{1:>20}".format("Bisulfite genome conversion start!",ctime()))
		now = datetime.now()
		genes = self.readGenome()

		c = 0
		hands = []
		slicenum = 0
		lengh = sumlength = 0
		for index,mtype in enumerate(DBFormatter.TYPES):
			#outfile = open(self.oname + "." + mtype + "_" + '{:04d}'.format(slicenum),"w")
			outfile = open(self.oname + "." + mtype,"w")
			hands.append(outfile)
		#moutfile = open(self.oname + ".MCT","w")
		#rmoutfile = open(self.oname + ".MGA","w")
		#uoutfile = open(self.oname + ".UCT","w")
		#ruoutfile = open(self.oname + ".UGA","w")

		for g in genes:
			id = g.id
			#print(g.id,len(g.seq))
			gle = len(g.seq)
			rawseq = str(g.seq).upper()
			tempseqs = ["","","",""]

			newseq = re.sub(r'[^AGCTN]','N',rawseq)
			tempseqs[0] = newseq.replace('CG','ZG')
			tempseqs[0] = tempseqs[0].replace('C','T')
			tempseqs[0] = tempseqs[0].replace('Z','C')
			tempseqs[1] = newseq.replace('CG','CZ')
			tempseqs[1] = tempseqs[1].replace('G','A')
			tempseqs[1] = tempseqs[1].replace('Z','G')

			tempseqs[2] = newseq.replace('C','T')
			tempseqs[3] = newseq.replace('G','A')

			#if sumlength > self.slicesize:
			#	for i in range(len(hands)):
			#		hands[i].close()
			#	slicenum += 1
			#	for index,mtype in enumerate(DBFormatter.TYPES):
			#		outfile = open(self.oname + "." + mtype + "_" + '{:04d}'.format(slicenum),"w")
			#		hands[index] = outfile
			#	sumlength = 0

			for index,type in enumerate(DBFormatter.TYPES):
				g.seq = Seq(tempseqs[index])
				g.id = type + "_" + id
				hands[index].write(g.format("fasta"))

			sumlength += gle


			#if c > 3 :break
		print("Convertion Done.")
		now1 = datetime.now()
		for i in hands:
			i.close()
		print("Cost time:" + str((now1-now).seconds) + " s")

	def indexdb(self):
		now = datetime.now()
		print("IndexDB start: " + ctime())
		for i,t in enumerate(DBFormatter.TYPES):
			print(str(i + 1) + "/" + str(len(DBFormatter.TYPES)) + " : " + self.oname + "." + t)
			#cmd = "makeblastdb -in " + self.oname + "." + t + " -dbtype nucl"
			cmd = "bowtie-build " + self.oname + "." + t + " --threads " + str(self.threads) + " " + self.oname + '.' + t
			#cmd = "hisat2-build " + self.oname + "." + t + " -p " + str(self.threads) + " " + self.oname + "." + t
			#print(cmd)
			res = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
			if res.returncode == 0:
				print("Success hisat2:",self.oname + "." + t)
			else:
				print("Error hisat2:",self.oname + "." + t)
			cmd1 = "samtools faidx " + self.oname + "." + t 
			res = subprocess.run(cmd1,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
			if res.returncode == 0:
				print("Success index:",self.oname + "." + t)
			else:
				print("Error index:",self.oname + "." + t)
			now1 = datetime.now()
			print("Cost time:" + str((now1-now).seconds) + " s")
		now1 = datetime.now()
		print("IndexDB Done" + ctime())
		print("All cost time:" + str((now1-now).seconds) + " s")

	def readGenome(self):
		print("Read and convert sequence...")
		for gene in SeqIO.parse(self.gfile,"fasta"):
			gene.description = ""
			yield gene

