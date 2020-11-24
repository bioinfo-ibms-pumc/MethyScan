#########################################################################
# File Name: searchdb.py
# > Author: CaoYinghao
# > Mail: caoyinghao@gmail.com 
# Created Time: Sun Nov  1 20:09:56 2020
#########################################################################
#! /usr/bin/python

import glob
import sys
import re
import os
from time import ctime
from Bio.SeqUtils import GC
from datetime import datetime
from Bio.Seq import Seq
import collections
import subprocess 
from utils import PrimerDesign
import matplotlib
import matplotlib.pyplot as plt

class Searcher(object):
    TYPES = ["MCT","MGA","UCT","UGA"]
    def __init__(self,inputfile,dbfile,outputfile,threads=4,searchtype = "ALL",mismatch=3,minlen=100,maxlen=500):
        self.ifile = inputfile
        self.dbfile = dbfile
        self.ofile = outputfile
        self.mismatch = mismatch
        self.threads = threads
        self.primers = collections.defaultdict(list)
        self.minlen = int(minlen)
        self.maxlen = int(maxlen)
        self.seqstr = ""
        self.searchtype = searchtype
        with open(self.ifile,"r") as infile:
            seqstr = ""
            for line in infile:
                if line.strip().startswith("#"):continue
                parts = line.strip().split("\t")
                self.primers[parts[0]].append(parts[1])
                self.primers[parts[0]].append(parts[2])
                seqstr += parts[0] + "\t" + parts[1] + "\t" + "I" * len(parts[1]) + "\t" + parts[2] + "\t" + "I" * len(parts[2]) + "\n"
            self.seqstr = seqstr.rstrip("\n")
        self.abb = {"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["A","C"],
                "B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","C","G","T"]}

    def checkdb(self):
        if int(self.mismatch) > 3:
            print("Error: mismatch should be no larger than 3")
            exit()
        #for i in Searcher.TYPES:
        #    if not os.path.exists(self.dbfile + "." + i):
        #        print("Error: " + self.dbfile + "." + i + " is not exist!")
        #        exit()
        return ""

    def alignment(self):
        self.checkdb()
        #self.domainjob()
        self.plotPrimers()

    def plotPrimers(self):
        matplotlib.use('AGG')
        plt.rcParams['font.family'] = u'Times New Roman'
        #fig = plt.figure(figsize=(6,6))
        strs = ""
        with open(self.ofile,"r") as infile:
            for line in infile:
                if line.strip().find("Amplicon detail") > -1:break
                strs += line
        primers = strs.split("-"*100)[:-1]
        matplotlib.use('AGG')
        fig = plt.figure(figsize=(6,6 + len(primers) / 10))
        ax = plt.subplot(111)
        c = 0
        colors = {"MCT":"steelblue","MGA":"red","UCT":"brown","UGA":"cyan"}
        sizes = {"6":1,"5":10,"4":50,"3":100,"2":150,"1":200,"0":300}
        ticknames = collections.defaultdict(list)
        for p in primers:
            names = re.search("PrimerName:\s+(.*)\n",p)
            res = p.split("TemplateResults:")[1]
            name = names.group(1).split("|")[-1]
            if len(name.split("_")) >3: 
                slices = name.split("_")
                name = slices[0] + "_" + slices[2] + "_" + slices[3]
            c += 1
            print(len(res.split("\n")))
            tempn = 0
            for line in res.split("\n"):
                tempn += 1
                if line == "":continue
                regex = re.search("Chr:\s(.*?)_.*Length:\s+(\S+).*Fmismatch:\s+(\S+).*Rmismatch:\s+(\S+)",line)
                mytype = regex.group(1)
                length = int(regex.group(2))
                mismatch = int(regex.group(3)) + int(regex.group(4))
                #print(mytype,length,mismatch)
                #print("hello",line)
                #print(mytype,length,mismatch)
                if int(mismatch) > self.mismatch * 2:continue
                if tempn > 50:break
                ticknames[name].append([length,str(mismatch),mytype])
        ax.set_xlim(self.minlen-50,self.maxlen+150)
        labels = []
        for c,name in enumerate(sorted(ticknames,key=lambda x:self.smallsort(x))):
            labels.append(re.sub("CT","SP",name))
            #labels.append(name)
            templen = -1
            for locs in sorted(ticknames[name]):
                length,mis,mytype = locs
                if templen == -1:
                    templen = length
                else:
                    if templen == length:
                        length += 2
                    else:
                        templen = length
                ax.scatter(y=c,x=length,marker="o",s=sizes[mis],edgecolor=colors[mytype],color="")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.text(self.minlen,len(ticknames) + 1,s= "Mismatches")
        for i in sizes:
            ax.scatter(self.minlen+200 + int(i) * 40,len(ticknames) +1.2,marker="o",s=sizes[i],edgecolor="black",color="")
            ax.text(self.minlen+200 + int(i) * 40-5,len(ticknames) +0.2 - 0.1,s=i)
        start = self.minlen
        for i,k in enumerate(colors):
            ax.hlines(xmin=start + 50 * i,xmax = start + 50 * i + 50,y=len(ticknames) + 2.5,color=colors[k])
            ax.text(x=start + 50 * i + 65,y=len(ticknames) + 2.5 -0.35,s=k)
            start += 80
        ax.set_xlabel("Amplicon length",size=15)
        plt.subplots_adjust(left=0.3)
        plt.savefig("primer_res.pdf")
        
    def smallsort(self,x):
        if x.find("_") > -1:
            xs = x.split("_")
            return [xs[1],xs[0],xs[-1]]
        else:
            return x

    def domainjob(self):
        print("{0:40}{1:>20}".format("Evaluation start", ctime()))
        titles = collections.defaultdict(dict)
        mutseqs = collections.defaultdict(dict)
        if self.searchtype == "ALL":
            for k,mytype in enumerate(Searcher.TYPES):
                print(str(k+1) + "/" + str(len(Searcher.TYPES)) + " " + mytype + " alignment start:" + ctime())
                temptitles,tempmutseqs = self.dojob(mytype)
                for t in temptitles:
                    for v in temptitles[t]:
                        titles[t][v] = temptitles[t][v]
                        mutseqs[t][v] = tempmutseqs[t][v]

        else:
            if self.searchtype not in Searcher.TYPES:
                print("Error: " + self.searchtype + "." + i + " is not defined!")
                exit()
            else:
                print("1/1 " + self.searchtype + " alignment start:" + ctime())
                titles,mutseqs = self.dojob(self.searchtype)
                #print(self.searchtype + " done")
                #c = 0
                #strs = ""
                #with open(self.ofile + "." + self.searchtype,"r") as infile:
                #    for line in infile:
                #        if line.strip().startswith("@"):continue
                #        strs += line
                #        c += 1
                #        #if c > 10:break
                #titles,mutseqs = self.parse(strs,self.searchtype)

        laststr = ""
        laststr += "R:[AG] Y:[CT] S:[GC] W:[AT] K:[GT] M:[AC] B:[CGT] D:[AGT] H:[ACT] V:[ACG] N:[ACGT]\n"
        laststr += """
MCT: Methylated genome with C converted to T except for CG sites.
UCT: Methylated genome with C converted to T.
MGA: Methylated genome with G converted to A except for CG sites.
UGA: Methylated genome with G converted to A.\n"""
        laststr += "#" * 100 + "\n" + "{0:40}{1:^20}{2:40}".format("#"*40,"Amplicon summary","#"*40) + "\n" + "#" * 100 + "\n"
        structs = collections.defaultdict(dict)
        for pp in sorted(titles):
            laststr += "PrimerName: " + pp + "\n"
            structs[pp]['TM'] = PrimerDesign.tmvalue(self.primers[pp][0])
            structs[pp]['GC'] = int(GC(Seq(self.primers[pp][0])) * 100)/100
            structs[pp]['SD'] = PrimerDesign.pcalcHomodimer(self.primers[pp][0])
            structs[pp]['SH'] = PrimerDesign.pcalcHairpin(self.primers[pp][0])
            structs[pp]['TM1'] = PrimerDesign.tmvalue(self.primers[pp][1])
            structs[pp]['GC1'] = int(GC(Seq(self.primers[pp][1])) * 100)/100
            structs[pp]['SD1'] = PrimerDesign.pcalcHomodimer(self.primers[pp][1])
            structs[pp]['SH1'] = PrimerDesign.pcalcHairpin(self.primers[pp][1])
            structs[pp]['HD'] = PrimerDesign.pcalcHeteDimer(self.primers[pp][0],self.primers[pp][1])
            laststr += " " + "{0:<30}{1:10}{2:10}".format(self.primers[pp][0],"TM:" + str(structs[pp]['TM']), "GC:"+str(structs[pp]['GC'])) + "\n"
            laststr += " " + "{0:<30}{1:10}{2:10}".format(self.primers[pp][1],"TM:" + str(structs[pp]['TM1']), "GC:"+str(structs[pp]['GC1'])) + "\n"
            if structs[pp]['SD'] != "" or structs[pp]['SH'] != "" or structs[pp]['SD1'] != "" or structs[pp]['SH1'] != "" or structs[pp]['HD'] != "":
                laststr += " " + "Dimer/Hairpin Found" + "\n\n"
            amps = titles[pp]
            if len(amps) == 0:
                laststr += "TemplateResults: None\n"
            else:
                laststr += "TemplateResults:\n"
            c = 0
            for j in sorted(amps,key = lambda x:[amps[x][0],amps[x][1],amps[x][2]]):
                c += 1
                laststr += "{0:15}{1}".format("Amplicon " + str(c),j) + "\n"
                if c > 1000:continue
            laststr += "-" * 100 + "\n\n"
        laststr += "\n" + "#" * 100 + "\n" + "{0:40}{1:^20}{2:40}".format("#"*40,"Amplicon detail","#"*40) + "\n" + "#" * 100 + "\n"
        for pp in sorted(titles):
            laststr += "PrimerName: " + pp + "\n"
            laststr += " " + "{0:<30}{1:10}{2:10}{3} {4}".format(self.primers[pp][0],"TM:" + str(structs[pp]['TM']), 
                    "GC:"+str(structs[pp]['GC']),structs[pp]['SD'],structs[pp]['SH']) + "\n"
            laststr += " " + "{0:<30}{1:10}{2:10}{3} {4}".format(self.primers[pp][1],"TM:" + str(structs[pp]['TM1']),
                    "GC:"+str(structs[pp]['GC1']),structs[pp]['SD1'],structs[pp]['SH1']) + "\n"
            laststr += " " + "{0:30}{1}".format("",structs[pp]['HD']) + "\n\n"
            amps = titles[pp]
            if len(amps) == 0:
                laststr += "TemplateResults: None\n"
            else:
                laststr += "TemplateResults:\n"
            c = 0
            for j in sorted(amps,key = lambda x:[amps[x][0],amps[x][1],amps[x][2]]):
                c += 1
                laststr += "{0:15}{1}".format("Amplicon " + str(c),j) + "\n"
                if c > 1000:continue
            laststr += "-" * 100 + "\n\n"
            laststr += "TemplateDetail:\n\n"
            c = 0
            for j in sorted(amps,key = lambda x:[amps[x][0],amps[x][1],amps[x][2]]):
                c += 1
                laststr += "{0:15}{1}".format("Amplicon " + str(c),j) + "\n"
                laststr += mutseqs[pp][j] + "\n"
                if c > 1000:continue
            laststr += "#" * 100 + "\n"
        if self.ofile != "":
            out = open(self.ofile,"w")
            out.write(laststr)
            out.close()
        else:
            print(laststr)


    def dojob(self,mytype):
        cmd = "echo \"" + self.seqstr + "\"| bowtie -y -r -e 90 -n 3 -l 5 -a " + self.dbfile + "." + mytype + " --sam-nohead -S -X 1000 -p " + str(self.threads) + " --12 -"
        #cmd = "echo \"" + self.seqstr + "\"| bowtie -y -r -e 90 -n 3 -l 5 -a " + self.dbfile + "." + mytype + " --sam-nohead -S mytest.sam -X 1000 -p " + str(self.threads) + " --12 -"
        #hisat2 -f -a --no-head -x " + self.dbfile + "." + i + " -p " + str(self.threads) + " -U " + self.ifile + " --very-sensitive"
        #print(cmd)
        res = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
        restitles = {}
        resmutseqs = {}
        if res.returncode == 0:
            #print(res)
            #out = open("cdkn." + i + ".sam","w")
            #out.write(res.stdout)
            #out.close()
            #exit()
            #print(res.stdout)
            #exit()
            restitles,resmutseqs = self.parse(res.stdout,mytype)
        else:
            print("Error alignment:",self.dbfile + "." + i,self.ifile)

        return restitles,resmutseqs

        pass

    def parse(self,samstr,mytype):
        print(" extract alignment start:" + ctime())
        sams = samstr.split("\n")
        dicts = collections.defaultdict(dict)
        for sam in sams:
            if sam == "":continue
            tokens = sam.split("\t")
            #print(sam)
            #exit()
            s = int(tokens[1])
            if not s & 4:
                misstag = tokens[-2].split(":")[-1]
                if int(misstag) > self.mismatch:continue
                isize = str(abs(int(tokens[8])))
                if abs(int(tokens[8])) > self.maxlen or abs(int(tokens[8])) < self.minlen:continue
                mp = MapLoc(loc = tokens[3],miss=misstag,isize=int(tokens[8]),length=len(tokens[9]))
                #print(isize,sam,abs(int(tokens[8])))
                if tokens[2] not in dicts[tokens[0]]:
                    dicts[tokens[0]][tokens[2]] = {}
                if isize not in dicts[tokens[0]][tokens[2]]:
                    dicts[tokens[0]][tokens[2]][isize] = []
                dicts[tokens[0]][tokens[2]][isize].append(mp)
        if len(dicts) > 0:
            print("Amplicons found for",len(dicts),"primer pairs!")
        else:
            print("No amplicons found!")
        amps = {} #collections.defaultdict(dict)
        for i in dicts:
            chrs = dicts[i]
            matchpairs = []
            for c in chrs:
                isizes = chrs[c]
                #print(isizes)
                for isize,lists in isizes.items():
                    #print(isize,lists)
                    plus = 0
                    pluspairs = {}
                    minuspairs = {}
                    for l in lists:
                        if l.size > 0:
                            pluspairs[l.loc] = l
                        else:
                            minuspairs[l.loc] = l
                    for loc in minuspairs:
                        loc1 = loc + minuspairs[loc].size + minuspairs[loc].length
                        if loc1 in pluspairs:
                            if mytype in ['MGA','UGA']:
                                amp = Amplicon(c,minuspairs[loc],pluspairs[loc1])
                            else:
                                amp = Amplicon(c,pluspairs[loc1],minuspairs[loc])
                            matchpairs.append(amp)
            amps[i]=matchpairs
        name = "loctime." + str(datetime.now().timestamp()) + ".bed"
        out = open(name,"w")
        fprimer = rprimer = ""
        cls = collections.defaultdict(dict)
        #print(amps)
        for amp in amps:
            #print(amp,self.primers[amp])
            #print(amps[amp])
            fprimer = self.primers[amp][0]
            rprimer = self.primers[amp][1]
            bedstr = ""
            for m in amps[amp]:
                if mytype in ['MGA','UGA']:
                    tempstr = str(m.chrom) + ":" + str(m.rloc.loc-1) + "-" + str(m.floc.loc + len(fprimer) - 1) + "(-)"
                    bedstr += str(m.chrom) + "\t" + str(m.rloc.loc-1) + "\t" + str(m.floc.loc + len(fprimer) - 1) + "\t\t\t-\n"
                    cls[amp][tempstr] = ""
                else:
                    tempstr = str(m.chrom) + ":" + str(m.floc.loc-1) + "-" + str(m.rloc.loc + len(rprimer) - 1)
                    bedstr += str(m.chrom) + "\t" + str(m.floc.loc-1) + "\t" + str(m.rloc.loc + len(rprimer) - 1) + "\t\t\t+\n"
                    cls[amp][tempstr] = ""
            #print("BEDSTR",bedstr)
            out.write(bedstr)
            #break
        out.close()
        cmd = "bedtools getfasta -fi " + self.dbfile.split(".")[0] + "." + mytype + " -bed " + name 
        if mytype in ["MGA","UGA"]:
            cmd += " -s"
        #print(cmd)
        #cmd = "ls "
        res = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8")
        #print(res.stdout.readline())
        #res.wait()
        restitles = collections.defaultdict(dict)
        resmutseqs = collections.defaultdict(dict)
        if res.returncode == 0:
            print("Success get sequences:")
            seqs = res.stdout.split(">")
            clsseq = {}
            for seq in seqs[1:]:
                title,myseq = seq.strip().split("\n")
                clsseq[title] = myseq
            #print(cls)
            #print(clsseq)
            for amp in cls:
                fprimer = self.primers[amp][0]
                rprimer = self.primers[amp][1]
                for tempstr in cls[amp]:
                    if tempstr in clsseq:
                        myseq = clsseq[tempstr]
                        #print(amp)
                        #print(fprimer)
                        #print(rprimer)
                        #print(tempstr)
                        #print(clsseq[tempstr])
                        #exit()
                        mychr,myloc = tempstr.split(":")
                        start,end = myloc.split("(")[0].split("-")
                        start = int(start)
                        end = int(end)
                        if mytype in ["MGA","UGA"]:
                            start,end = end,start
                        newfp = []
                        newrp = []
                        lenf = len(fprimer)
                        lenr = len(rprimer)
                        lenall = len(myseq)
                        if mytype in ['MCT','UCT']:
                            primer1end = start + lenf
                            primer2end = end - lenr + 1
                            start += 1
                            ftemp = myseq[:len(fprimer)]
                            rtemp = Seq(myseq[-len(rprimer):]).reverse_complement()
                        else:
                            primer1end = start - lenf + 1
                            primer2end = end + lenr
                            end += 1
                            ftemp = myseq[:len(fprimer)]
                            rtemp = Seq(myseq[-len(rprimer):]).reverse_complement()
                        fmis = 0
                        rmis = 0
                        for i in range(lenf):
                            if ftemp[i] == fprimer[i]:
                                newfp.append("|")
                            else:
                                if fprimer[i] in self.abb:
                                    if ftemp[i] in self.abb[fprimer[i]]:
                                        newfp.append("~")
                                    else:
                                        newfp.append(" ")
                                        fmis += 1
                                else:
                                    newfp.append(" ")
                                    fmis += 1
                        for i in range(lenr):
                            if rtemp[i] == rprimer[i]:
                                newrp.append("|")
                            else:
                                if rprimer[i] in self.abb:
                                    if rtemp[i] in self.abb[rprimer[i]]:
                                        newrp.append("~")
                                    else:
                                        newrp.append(" ")
                                        rmis += 1
                                else:
                                    newrp.append(" ")
                                    rmis += 1
                        
                        #mutseq = "".join(newfp) + "." * (lenall - lenf - lenr) + "".join(newrp)
                        #print(amp)
                        if fmis + rmis > self.mismatch * 2:continue
                        alen = end - start + 1 if end > start else start - end + 1

                        muttitle = "{0:20}{1:20}{2:20}{3:15}{4:15}".format("Chr: " + mychr,"Location: " + str(start), "Length: " + str(alen),"Fmismatch: " + str(fmis), "Rmismatch: " + str(rmis))
                        restitles[amp][muttitle] = [fmis + rmis,alen,mytype]
                        #print(muttitle)
                        mutseq = ""
                        mutseq += "{0:<15}: {1:>10} 5'- {2} -3' {3:<2d}".format("Forward primer","1",fprimer,lenf) + "\n"
                        mutseq += "{0:<32}{1}".format("","".join(newfp)) + "\n"
                        mutseq += "{0:<15}: {1:>10} 5'- {2} -3' {3:<2d}".format("Template",start,ftemp,primer1end) + "\n\n"
                        mutseq += "{0:<15}: {1:>10} 5'- {2} -3' {3:<2d}".format("Reverse primer","1",rprimer,lenr) + "\n"
                        mutseq += "{0:<32}{1}".format("","".join(newrp)) + "\n"
                        mutseq += "{0:<15}: {1:>10} 5'- {2} -3' {3:<2d}".format("RC_Template",end,rtemp,primer2end) + "\n\n"
                        mutseq += "{0:<15}: {1:>10} 5'- {2:} -3' {3}".format(mytype + "_template ",start,myseq,end) + "\n"
                        tm = PrimerDesign.tmvalue(myseq)
                        gc = int(GC(Seq(myseq)) * 100)/100
                        mutseq += "{0:<15} TM:{1:<5} GC:{2:<5} Length:{3:<5}".format(mytype + "_template ",tm,gc,alen) + "\n"
                        #print("{0:<10}:{1:300}".format("Primers",mutseq))
                        #print(mutseq)
                        #break
                        resmutseqs[amp][muttitle] = mutseq
        else:
            print("Error get sequences:")
        print(" extract alignment done:" + ctime())
        if os.path.exists(name):
            os.remove(name)
        return restitles,resmutseqs

            #print(cmd)

class Amplicon(object):
    def __init__(self,chrom,loc1,loc2):
        self.chrom = chrom
        self.floc = loc1
        self.rloc = loc2
        self.size = abs(loc1.size)
        self.miss = loc1.miss + loc2.miss
        pass

    def setforwardprimer(self,primer):
        self.fprimer = primer

    def __repr__(self):
        return "Chr: " + self.chrom + " Size: " + str(self.size) + " mismatch: " + str(self.miss) + " " + str(self.floc) + "\t" + str(self.rloc)

class MapLoc(object):
    def __init__(self,loc,miss,isize,length):
        self.loc = int(loc)
        self.miss = int(miss)
        self.size = isize
        self.length = length
        pass

    def __repr__(self):
        return " loc: " + str(self.loc) + " mismatch: " + str(self.miss) + " size: " + str(self.size) + " length: " + str(self.length)

    #def __str__(self):
    #    return "chr: " + self.chrom + " loc: " + str(self.loc) + " mismatch: " + self.miss + " size: " + str(self.size)
