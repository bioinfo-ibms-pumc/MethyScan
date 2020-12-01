#########################################################################
# File Name: utils.py
# > Author: CaoYinghao
# > Mail: caoyinghao@gmail.com 
# Created Time: Sat Oct 31 10:38:04 2020
#########################################################################
#! /usr/bin/python

import glob
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import collections
import primer3
import matplotlib
import matplotlib.pyplot as plt

class PrimerDesign(object):
    PRIMERTYPE = ("BSP","USP","MSP","ALL")
    def __init__(self,seqfile,threshold,type="BSP",maxnum = 3):
        self.seqfile = seqfile
        self.threshold = threshold
        self.primertype = type
        self.maxnum = maxnum

    def designprimers(self,outfile = ""):
        #if outfile == "":
        pmx = self.threshold.primermaxlen
        pmi = self.threshold.primerminlen
        results = ""
        pseq = ""
        for gene in SeqIO.parse(self.seqfile,"fasta"):
            seq = gene.seq
            if self.threshold.regionend == -1:
                self.threshold.regionend = len(seq)
            if self.threshold.regionstart >= len(seq):
                print("Error: Region start is larger than the sequence " + seq.id)
                exit()
            if self.threshold.regionstart >= self.threshold.regionend:
                print("Error: Region start is larger than region end")
                exit()
            if self.threshold.regionend - self.threshold.regionstart + 1 <= 100:
                print("Error: Amplicon length must be larger than 100 bp")
            ### nest
            results += self.threshold.print() 

            usedstr = str(seq[self.threshold.regionstart-1:self.threshold.regionend-1]).upper()
            nums = []
            nums1 = []
            for i in range(0,len(usedstr)-100,1):
                temps = usedstr[i:i + 100]
                dicts = collections.defaultdict(int)
                for j in temps:
                    dicts[j] += 1
                gc = GC(Seq(temps))/ 100
                cpgs = len(temps.split("CG"))-1
                oe = 0
                if dicts['C'] != 0 and dicts['G'] != 0:
                    oe = cpgs * len(temps) / dicts["C"] / dicts["G"]
                if oe > 0.6 and gc > 0.5:
                    nums.append(i)
                nums1.append(gc)
                #if i > 1500:break
            #print(nums)
            segments = []
            start = -1
            end = -1
            regions = []
            regionlocs = []
            ys = []
            for i in nums:
                if start == -1:
                    start = i
                    end = start
                    continue
                else:
                    if i - end < 30:
                        end = i
                    else:
                        #print(start,end)
                        if end - start > 100:
                            locs = []
                            ys = []
                            regionlocs.append((start+50,end+50))
                            for k in range(start+50,end+50):
                                locs.append(k)
                                ys.append(nums1[k-50])
                            regions.append([locs,ys])
                        start = i
                        end = i
            if end - start > 100:
                locs = []
                ys = []
                regionlocs.append((start+50,end+50))
                for k in range(start+50,end+50):
                    locs.append(k)
                    ys.append(nums1[k-50])
                regions.append([locs,ys])
            #print(regionlocs)

            matplotlib.use('AGG')
            plt.rcParams['font.family'] = u'Times New Roman'
            fig = plt.figure(figsize=(15,5))
            plt.plot(range(50,len(nums1)+50),nums1,"--",color="black")
            #print(len(ys))
            for region in regions:
                plt.stackplot(region[0],region[1],color="cyan",alpha=0.5)
            plt.ylim(0,1)
            lens = len(usedstr) / 7
            lwd = 0.5
            plt.hlines(xmin=lens,xmax=lens * 1.25,y=0.9,color="red",linewidth=lwd*10)
            plt.hlines(xmin=lens*1.75,xmax=lens * 2,y=0.9,color="red",linewidth=lwd*10)
            plt.hlines(xmin=lens*1.25,xmax=lens * 1.75,y=0.9,color="gray",linewidth=lwd*5,linestyle="dotted")
            plt.text(x=lens*2.25,y=0.89,s="M primers")
            plt.hlines(xmin=lens*3,xmax=lens * 3.25,y=0.9,color="steelblue",linewidth=lwd*10)
            plt.hlines(xmin=lens*3.75,xmax=lens * 4,y=0.9,color="steelblue",linewidth=lwd*10)
            plt.hlines(xmin=lens*3.25,xmax=lens * 3.75,y=0.9,color="gray",linewidth=lwd*5,linestyle="dotted")
            plt.text(x=lens*4.25,y=0.89,s="U primers")
            plt.hlines(xmin=lens*5,xmax=lens * 5.25,y=0.9,color="black",linewidth=lwd*10)
            plt.hlines(xmin=lens*5.75,xmax=lens * 6,y=0.9,color="black",linewidth=lwd*10)
            plt.hlines(xmin=lens*5.25,xmax=lens * 5.75,y=0.9,color="gray",linewidth=lwd*5,linestyle="dotted")
            plt.text(x=lens*6.25,y=0.89,s="Nested primers")
            #exit()
            #for i in regionlocs:
            #cpglocs = regionlocs[-1]

            if len(regionlocs) == 0:
                print("No CpG islands found in " + gene.id)
                regionlocs = [[1,len(usedstr)]]
            else:
                print(len(regionlocs),"CpG islands found in " + gene.id)
            for index,cpgloc in enumerate(regionlocs):
                #if index != 1:continue
                res,aseq,mseq,useq = self.calcu(usedstr,cpgloc)
                h = 0.05 * len(res) * (index % 2)
                results += "\n" + "#" * 20 + "{0:^55}".format(gene.id + " CpG islands " + str(index + 1) + " location: " + str(cpgloc[0]) + "-" + str(cpgloc[1])) + "#" * 20 + "\n\n"
                results += "{0:<30}: ".format("Reference sequence") + aseq + "\n"
                results += "{0:<30}: ".format("Bisulfite modified sequence") + mseq + "\n"
                vseq = re.sub("CG","++",mseq)
                vseq = re.sub("\w",".",vseq)
                results += "{0:<18}{1:<3}{2:<9}: ".format("Symbols (CpG num: ",len(mseq.split("CG"))-1,")") + vseq + "\n"
                temppseq = ""
                for index1,r in enumerate(res):
                    cpgs,mfp,mrp,ufp,urp,nfp,nrp,allcpgs = r
                    wseq = re.sub(mfp[2],">" * len(mfp[2]),mseq)
                    wseq = re.sub(str(Seq(mrp[2]).reverse_complement()),"<" * len(mrp[2]),wseq)
                    wseq = re.sub("\w",".",wseq)
                    results += "{0:<25}{1:<2}{2:<3}: ".format("M primers " + str(index1 + 1) + " (CpG num: ",cpgs,")")  + wseq + "\n"
                    wseq = re.sub(ufp[2],"|" * len(ufp[2]),useq)
                    wseq = re.sub(str(Seq(urp[2]).reverse_complement()),"|" * len(urp[2]),wseq)
                    wseq = re.sub("\w",".",wseq)
                    results += "{0:<30}: ".format("U primers " + str(index1 + 1)) + wseq + "\n"
                    plt.hlines(xmin=mfp[1],xmax=mfp[1] + mfp[0],y=h + 0.034,color="red",linewidth=lwd*4)
                    plt.hlines(xmin=mrp[1],xmax=mrp[1] + mrp[0],y=h + 0.034,color="red",linewidth=lwd*4)
                    plt.hlines(xmin=mfp[1] + mfp[0],xmax=mrp[1],y=h + 0.034,color="gray",linewidth=lwd*2,linestyle="dotted")
                    plt.hlines(xmin=ufp[1],xmax=ufp[1] + ufp[0],y=h + 0.048,color="steelblue",linewidth=lwd*4)
                    plt.hlines(xmin=urp[1],xmax=urp[1] + urp[0],y=h + 0.048,color="steelblue",linewidth=lwd*4)
                    plt.hlines(xmin=ufp[1] + ufp[0],xmax=urp[1],y=h + 0.048,color="gray",linewidth=lwd,linestyle="dotted")
                    plt.hlines(xmin=nfp[1],xmax=nfp[1] + nfp[0],y=h + 0.02,color="black",linewidth=lwd*4)
                    plt.hlines(xmin=nrp[1],xmax=nrp[1] + nrp[0],y=h + 0.02,color="black",linewidth=lwd*4)
                    plt.hlines(xmin=nfp[1] + nfp[0],xmax=nrp[1],y=h + 0.02,color="gray",linewidth=lwd*2,linestyle="dotted")
                    h += 0.05
                    titlemf = "p" + str(index1 + 1) + "_" + str(mfp[1])+ "_CpG-" + str(cpgloc[0]) + "-" + str(cpgloc[1]) + "_M_F"
                    titlemr = "p" + str(index1 + 1) + "_" + str(mrp[1])+ "_CpG-" + str(cpgloc[0]) + "-" + str(cpgloc[1]) + "_M_R"
                    titleuf = "p" + str(index1 + 1) + "_" + str(ufp[1])+ "_CpG-" + str(cpgloc[0]) + "-" + str(cpgloc[1]) + "_U_F"
                    titleur = "p" + str(index1 + 1) + "_" + str(urp[1])+ "_CpG-" + str(cpgloc[0]) + "-" + str(cpgloc[1]) + "_U_R"
                    titlenf = "p" + str(index1 + 1) + "_" + str(nfp[1])+ "_CpG-" + str(cpgloc[0]) + "-" + str(cpgloc[1]) + "_Nested_F"
                    titlenr = "p" + str(index1 + 1) + "_" + str(nrp[1])+ "_CpG-" + str(cpgloc[0]) + "-" + str(cpgloc[1]) + "_Nested_R"
                    pseq += gene.id + "|" + titlemf[:-2] + "\t" + mfp[2] +"\t" + mrp[2] + "\n"
                    pseq += gene.id + "|" + titleuf[:-2] + "\t" + ufp[2] +"\t" + urp[2] + "\n"
                    pseq += gene.id + "|" + titlenf[:-2] + "\t" + nfp[2] +"\t" + nrp[2] + "\n"
                    temppseq += "{0:>40}{1:<30}{2:<10}{3:<10}{4:5}{5:5}{6:<10}\n".format(titlemf + "orward ","Sequence: " + mfp[2], " Length: " + str(mfp[0]),
                            " Location: " + str(mfp[1])," TM: " + str(mfp[-2])," GC: " + str(mfp[-1])," CpG Num: " + str(mfp[3]))
                    temppseq += "{0:>40}{1:<30}{2:<10}{3:<10}{4:5}{5:5}{6:<10}\n".format(titlemr + "everse ","Sequence: " + mrp[2], " Length: " + str(mrp[0]),
                            " Location: " + str(mrp[1])," TM: " + str(mrp[-2])," GC: " + str(mrp[-1])," CpG Num: " + str(mrp[3]))
                    temppseq += "{0:>40}{1:<30}{2:<10}{3:<10}{4:5}{5:5}{6:<10}\n".format(titleuf + "orward ","Sequence: " + ufp[2], " Length: " + str(ufp[0]),
                            " Location: " + str(ufp[1])," TM: " + str(ufp[-2])," GC: " + str(ufp[-1])," CpG Num: " + str(ufp[3]))
                    temppseq += "{0:>40}{1:<30}{2:<10}{3:<10}{4:5}{5:5}{6:<10}\n".format(titleur + "everse ","Sequence: " + urp[2], " Length: " + str(urp[0]),
                            " Location: " + str(urp[1])," TM: " + str(urp[-2])," GC: " + str(urp[-1])," CpG Num: " + str(urp[3]))
                    temppseq += "{0:>40}{1:<30}{2:<10}{3:<10}{4:5}{5:5}{6:<10}\n".format(titlenf + "orward ","Sequence: " + nfp[2], " Length: " + str(nfp[0]),
                            " Location: " + str(nfp[1])," TM: " + str(nfp[-2])," GC: " + str(nfp[-1])," CpG Num: " + str(nfp[3]))
                    temppseq += "{0:>40}{1:<30}{2:<10}{3:<10}{4:5}{5:5}{6:<10}\n".format(titlenr + "everse ","Sequence: " + nrp[2], " Length: " + str(nrp[0]),
                            " Location: " + str(nrp[1])," TM: " + str(nrp[-2])," GC: " + str(nrp[-1])," CpG Num: " + str(nrp[3]))
                    ampcon = usedstr[mfp[1]:mrp[1]]
                    temppseq += "{0:>49} {1}".format("M Amplicon Length:" , str(mrp[1] - mfp[1] )) + "\n"
                    temppseq += "{0:>49} {1}".format("U Amplicon Length:" , str(urp[1] - ufp[1])) + "\n"
                    temppseq += "{0:>49} {1}".format("Nested Amplicon Length:" , str(nrp[1] - nfp[1])) + "\n"
                    temppseq += "\n"
                    #break
                if len(res) == 0:
                    results += "{0:<40}".format("No primer in this region!") + "\n"
                else:
                    results += temppseq + "\n"
                    pass
                results += "\n" + "#" * 95 + "\n\n"
                #break
            plt.xlabel("Sequence")
            plt.ylabel("GC Percentage")
            plt.yticks([0,0.2,0.4,0.6,0.8,1],["0","20","40","60","80","100"])

            if outfile is not None:
                plt.savefig(outfile + "_" + gene.id + ".pdf")
            else:
                plt.savefig(gene.id + ".pdf")
        if outfile is not None:
            out = open(outfile + ".fa","w")
            out.write(pseq)
            out.close()
            out = open(outfile + ".detail.txt","w")
            out.write(results)
            out.close()
        else:
            print("\nPrimer results\n")
            print(pseq)
            print(results)

            #for nf in range(self.threshold.regionstart-1,self.threshold.regionend - 1):
            #    for fmx in range(self.threshold)

    @staticmethod
    def tmvalue(sequence):
        #tm = mt.Tm_NN(Seq(sequence),Na=50,Mg=0,dnac1=250,dnac2=250)
        if sequence.find("R") > -1:
            sequence = re.sub("R","G",sequence)
        if sequence.find("Y") > -1:
            sequence = re.sub("Y","C",sequence)
        tm = primer3.calcTm(sequence,tm_method="breslauer",salt_corrections_method="schildkraut")
        return int(tm * 100)/100

    @staticmethod
    def pcalcHomodimer(seq):
        res = primer3.calcHomodimer(seq)
        dg = "{0:.2f}".format(res.dg)
        tm = "{0:.2f}".format(res.tm)
        if res.structure_found == True and float(dg) < -1:
            return "{0}{1}".format("Self_Dimer:", "+:tm:" + tm + ";deltaG:" + dg)
        else:
            return ""

    @staticmethod
    def pcalcHeteDimer(seq1,seq2):
        res = primer3.calcHeterodimer(seq1,seq2)
        dg = "{0:.2f}".format(res.dg)
        tm = "{0:.2f}".format(res.tm)
        if res.structure_found == True and float(dg) < -1:
            return "{0}{1}".format("HeteDimer:", "+:tm:" + tm + ";deltaG:" + dg)
        else:
            return ""

    @staticmethod
    def pcalcHairpin(seq):
        res = primer3.calcHairpin(seq)
        dg = "{0:.2f}".format(res.dg)
        tm = "{0:.2f}".format(res.tm)
        if res.structure_found == True and float(dg) < -1:
            return "{0}{1}".format("Self_Hairpin:", "+:tm:" + tm + ";deltaG:" + dg)
        else:
            return ""

    def calcuDH(self,m1,m2,u1,u2,n1,n2):
        count = 0
        if PrimerDesign.pcalcHairpin(m1) != "": count += 1
        if PrimerDesign.pcalcHairpin(m2) != "": count += 1
        if PrimerDesign.pcalcHairpin(u1) != "": count += 1
        if PrimerDesign.pcalcHairpin(u2) != "": count += 1
        if PrimerDesign.pcalcHairpin(n1) != "": count += 1
        if PrimerDesign.pcalcHairpin(n2) != "": count += 1

        if PrimerDesign.pcalcHomodimer(m1) != "": count += 1
        if PrimerDesign.pcalcHomodimer(m2) != "": count += 1
        if PrimerDesign.pcalcHomodimer(u1) != "": count += 1
        if PrimerDesign.pcalcHomodimer(u2) != "": count += 1
        if PrimerDesign.pcalcHomodimer(n1) != "": count += 1
        if PrimerDesign.pcalcHomodimer(n2) != "": count += 1

        if PrimerDesign.pcalcHeteDimer(m1,m2) != "": count += 1
        if PrimerDesign.pcalcHeteDimer(u1,u2) != "": count += 1
        if PrimerDesign.pcalcHeteDimer(n1,n2) != "": count += 1
        return count

    def calcu(self,usedstr,locs):
        ampseq = usedstr[locs[0]:locs[1]]
        mampseq = re.sub("CG","MG",ampseq)
        mampseq = re.sub("C","T",mampseq)
        mampseq = re.sub("M","C",mampseq)

        uampseq = re.sub("C","T",ampseq)
        #print(ampseq)
        #print(mampseq)
        #print(uampseq)
        rmampseq = str(Seq(mampseq).reverse_complement())
        uampseq = re.sub("C","T",mampseq)
        ruampseq = str(Seq(uampseq).reverse_complement())
        fps = []
        rps = []
        #print(locs)
        for i in range(self.threshold.primerminlen,self.threshold.primermaxlen+1):
            for j in range(0,len(mampseq)-i,1):
                s = mampseq[j:j+i]
                #if j + locs[0] == 695:
                #    print(s,j + locs[0])
                if s[-3:].find("CG") > -1:
                    fps.append((i,j + locs[0],s,len(s.split("CG"))-1,PrimerDesign.tmvalue(s),int(GC(Seq(s)) * 100)/100))
        for i in range(self.threshold.primerminlen,self.threshold.primermaxlen+1):
            for j in range(0,len(rmampseq)-i,1):
                s = rmampseq[j:j+i]
                #if len(rmampseq) + locs[0] - j == 814:
                #    print(i,j,s)
                if s[-3:].find("CG") > -1:
                    rps.append((i,len(rmampseq) + locs[0] - j ,s,len(s.split("CG"))-1,PrimerDesign.tmvalue(s),int(GC(Seq(s))*100)/100))

        candidates = {}
        for fp in fps:
            if fp[2].find("T"*self.threshold.polyT) > -1:continue
            if fp[-2] > self.threshold.maxtm or fp[-2] < self.threshold.mintm:continue

            for rp in rps:
                if rp[2].find("T"*self.threshold.polyT) > -1:continue
                if rp[-2] > self.threshold.maxtm or rp[-2] < self.threshold.mintm:continue
                if rp[1] - fp[1] < self.threshold.ampminlen:continue
                if rp[1] - fp[1] > self.threshold.ampmaxlen:continue

                cgs = rp[-3] + fp[-3]
                tmdiff = abs(fp[-2] - rp[-2])
                if tmdiff > self.threshold.tmdiff:continue
                #print(fp,rp,tmdiff)
                candidates[str(fp[0]) + "_" + str(fp[1]) + "_" + str(fp[2]) + "|" + str(rp[0]) +"_" + str(rp[1]) + "_" + str(rp[2])] = (cgs,tmdiff,rp[1] - fp[1],abs(fp[0] - rp[0]),fp,rp)
        ucandidates = self.refineprimer(candidates,uampseq,ruampseq,locs[0])
        nestcandidates = self.refinenest(ucandidates,usedstr,locs)
        #print(nestcandidates['22_658_TGTTTTTTTTTTAGGAAATAAC|25_916_TCATTCATTATACTACAATCTAACA'])
        #print(nestcandidates)
        print("M,U,Nest:",len(candidates),len(ucandidates),len(nestcandidates))
        names = []
        dimer_hairpin = {}
        for i in ucandidates:
            #if i != "20_662_TTTTTTTTAGGAAATAACGC|20_809_ATTATTAACGTTACCCTCGC":continue
            #print(i)
            if i in nestcandidates:
                dimer_hairpin[i] = self.calcuDH(candidates[i][-2][2],candidates[i][-1][2],
                        ucandidates[i][-2][2],ucandidates[i][-1][2],
                        nestcandidates[i][-3][2],nestcandidates[i][-2][2])
                names.append(i)
        print("Primers predicted in",len(names),"candidates in CpG island",locs)

        res = []
        c = 0
        for name in sorted(names,key=lambda x:self.sortfunc(x,candidates,ucandidates,nestcandidates,dimer_hairpin)):
            cpgs = candidates[name][0]
            mfp = candidates[name][-2]
            mrp = candidates[name][-1]
            ufp = ucandidates[name][-2]
            urp = ucandidates[name][-1]
            nfp = nestcandidates[name][-3]
            nrp = nestcandidates[name][-2]
            allcpgs = nestcandidates[name][-1]
            res.append([cpgs,mfp,mrp,ufp,urp,nfp,nrp,allcpgs])
            c += 1
            if c >= self.maxnum:break
            #print(name)
            #print(candidates[name])
            #print(ucandidates[name])
            #print(nestcandidates[name])
        return res,ampseq,mampseq,uampseq

    
    def sortfunc(self,x,candi,ucandi,ncandi,dimer_hairpin):
        cpgs = -candi[x][0]
        dh = dimer_hairpin[x]
        tmdiff = candi[x][1]
        tmdiff1 = max(candi[x][-2][-2],candi[x][-1][-2],ucandi[x][-2][-2],ucandi[x][-1][-2])-min(candi[x][-2][-2],candi[x][-1][-2],ucandi[x][-2][-2],ucandi[x][-1][-2])
        acpgs = ncandi[x][-1]
        return cpgs,dh,tmdiff,tmdiff1,acpgs


    def refinenest(self,candidates,usedstr,locs):
        mampseq = re.sub("CG","MG",usedstr)
        mampseq = re.sub("C","T",mampseq)
        mampseq = re.sub("M","C",mampseq)
        rmampseq = str(Seq(mampseq).reverse_complement())
        nestcandidates = {}
        #ampseq = usedstr[locs[0]:locs[1]]
        for v in sorted(candidates):
            ufp = candidates[v][0]
            urp = candidates[v][1]

            if ufp[1] - 25 < 0:
                print("Nest primer can't be predict for the location " + str(ufp[1]))
                exit()
            if urp[1] + 25 > len(usedstr):
                print("Nest primer can't be predict for the location " + str(urp[1]))
                exit()
            outside = 100
            start = ufp[1] - outside if ufp[1] - outside > 0 else 0
            length = len(usedstr)
            end = length - (urp[1] + outside - self.threshold.primerminlen) if urp[1] + outside < length else length - self.threshold.primerminlen
            #print(ufp[1],urp[1],locs)
            #print(start,end)
            #exit()
            fps = []
            rps = []
            for i in range(self.threshold.primerminlen,self.threshold.primermaxlen+1):
                for j in (start,ufp[1] - 25):
                    s = mampseq[j:j+i]
                    if s.find("T"*self.threshold.polyT) > -1:continue
                    #if j + locs[0] == 695:
                    #    print(s,j + locs[0])
                    temptm = PrimerDesign.tmvalue(s)
                    if temptm > self.threshold.maxtm or temptm < self.threshold.mintm:continue
                    if s.find("CG") == -1:
                        fps.append((i,j,s,len(s.split("CG"))-1,temptm,int(GC(Seq(s))*100)/100))
            for i in range(self.threshold.primerminlen,self.threshold.primermaxlen+1):
                for j in (length - (urp[1] + 25),end):
                    #s = rmampseq[j:j+i]
                    rs = rmampseq[j:j+i]
                    if rs.find("T"*self.threshold.polyT) > -1:continue
                    rtemptm = PrimerDesign.tmvalue(rs)
                    if rtemptm > self.threshold.maxtm or rtemptm < self.threshold.mintm:continue
                    #if len(rmampseq) + locs[0] - j == 814:
                    #    print(i,j,s)
                    if rs.find("CG") == -1:
                        rps.append((i,length - j,rs,len(rs.split("CG"))-1,rtemptm,int(GC(Seq(rs))*100)/100))
            tempcands = {}
            ###refind
            if len(fps) == 0:
                for i in range(self.threshold.primerminlen,self.threshold.primermaxlen+1):
                    for j in (start,ufp[1] - 25):
                        s = mampseq[j:j+i]
                        if s.find("T"*self.threshold.polyT) > -1:continue
                        temptm = PrimerDesign.tmvalue(s)
                        if temptm > self.threshold.maxtm or temptm < self.threshold.mintm:continue
                        if len(s.split("CG")) == 2:
                            s = re.sub("CG","YG",s)
                        if s.find("CG") == -1:
                            fps.append((i,j,s,1,temptm,int(GC(Seq(s))*100)/100))
            ###refind
            if len(rps) == 0:
                for i in range(self.threshold.primerminlen,self.threshold.primermaxlen+1):
                    for j in (length - (urp[1] + 25),end):
                        #s = rmampseq[j:j+i]
                        rs = rmampseq[j:j+i]
                        if rs.find("T"*self.threshold.polyT) > -1:continue
                        rtemptm = PrimerDesign.tmvalue(rs)
                        if rtemptm > self.threshold.maxtm or rtemptm < self.threshold.mintm:continue
                        if len(rs.split("CG")) == 2:
                            rs = re.sub("CG","CR",rs)
                        #if len(rmampseq) + locs[0] - j == 814:
                        #    print(i,j,s)
                        #print(rs)
                        if rs.find("CG") == -1:
                            rps.append((i,length - j,rs,1,rtemptm,int(GC(Seq(rs))*100)/100))

            for fp in fps:

                for rp in rps:
                    if rp[1] - fp[1] < 100:continue

                    cgs = len(usedstr[fp[1]:rp[1]].split("GC")) - 1
                    tmdiff = abs(fp[-2] - rp[-2])
                    if tmdiff > self.threshold.tmdiff:continue
                    tempcands[str(fp[0]) + "_" + str(fp[1]) + "_" + str(fp[2]) + "|" + str(rp[0]) +"_" + str(rp[1]) + "_" + str(rp[2])] = (cgs,tmdiff,rp[1] - fp[1],abs(fp[0] - rp[0]),fp,rp)
            for cand in sorted(tempcands,key = lambda x:(-tempcands[x][0],tempcands[x][1],tempcands[x][2],tempcands[x][3])):
                #print(tempcands[cand][-2],tempcands[cand][-1],locs)
                #exit()
                nestcandidates[v] = (tempcands[cand][-2],tempcands[cand][-1],tempcands[cand][0])
                #print(i)
                break
        return nestcandidates


    def refineprimer(self,candidates,ampseq,rampseq,shiftloc):
        num = 0
        ucandidates = {}
        for i in sorted(candidates):#,key=lambda x:(-candidates[x][0],candidates[x][1],candidates[x][2],candidates[x][3])):
            fp = candidates[i][-2]
            rp = candidates[i][-1]
            #if fp[1] != 696:
            #    continue
            #print(i,candidates[i])
            s = rs = loc = loc1 = tm = tm1= -1
            for j in range(0,self.threshold.primermaxlen +1- fp[0]):
                for k in range(0,self.threshold.primermaxlen +1 - fp[0]):
                    #print(j,k,fp,fp[1],len(uampseq),fp[1]-locs[0],fp[1] + k - locs[0])
                    if fp[1] - shiftloc - j > 0 and fp[1] - shiftloc+ k < len(ampseq):
                        s = ampseq[fp[1]-shiftloc -j :fp[1] -shiftloc + fp[0] + k]
                        if s.find("T"*self.threshold.polyT) > -1:continue
                        temptm = PrimerDesign.tmvalue(s)
                        if temptm > self.threshold.maxtm or temptm < self.threshold.mintm:continue
                        #print(s,temptm)
                        if abs(temptm  - fp[-2]) < self.threshold.tmdiff:
                            loc = fp[1] - j
                            tm = temptm
                            break
                if tm != -1:
                    #print("UF",loc,tm,loclen,s)
                    break
            for k in range(0,self.threshold.primermaxlen +1- rp[0]):
                for j in range(0,self.threshold.primermaxlen +1- rp[0]):
                    #print(rp[1], shiftloc,  j, len(rampseq),rp[0])
                    if rp[1] -shiftloc  + j < len(rampseq) and rp[1] -shiftloc + rp[0] + k < len(rampseq):
                        rs = rampseq[len(rampseq) - (rp[1] - shiftloc) - j :len(rampseq) - (rp[1] - shiftloc) + rp[0] + k]
                        #print(s,rs)
                        if rs.find("T"*self.threshold.polyT) > -1:continue
                        temptm = PrimerDesign.tmvalue(rs)
                        if temptm > self.threshold.maxtm or temptm < self.threshold.mintm:continue
                        if abs(temptm  - rp[-2]) < self.threshold.tmdiff:
                            loc1 = rp[1] - j
                            tm1 = temptm
                            ##print("hello",loc1,tm1,loclen1,rs)
                            break
                if tm1 != -1:
                    #print("UR",loc1,tm1,loclen1,rs)
                    break
            #print(s,rs)
            if tm != -1 and tm1 != -1:
                num += 1
                ufp = (len(s),loc,s,len(s.split("CG"))-1,tm,int(GC(Seq(s))*100)/100)
                urp = (len(rs),loc1,rs,len(s.split("CG"))-1,tm1,int(GC(Seq(rs))*100)/100)
                ucandidates[i] = (ufp,urp)
                #print("Welldone", s,loc,tm,rs,loc1,tm1,fp,rp,candidates[i][:4])
                #print(self.tmvalue(fp[2]))
                #print(primer3.calcTm("attttttaggtttcgtttcggc"))
                #print(primer3.calcTm("ttagttattttttaggttttgttttggt"))
                #break
        return ucandidates
                        
                    

class Threshold(object):
    def __init__(self,maxlen=25,minlen=20,ampmaxlen=200,ampminlen=100,maxtm=60.0,mintm=50.0,tmdiff=5.0,regionstart=1,regionend=-1):
        self.primermaxlen = int(maxlen)
        self.primerminlen = int(minlen)
        self.maxtm = float(maxtm)
        self.mintm = float(mintm)
        self.ampmaxlen = ampmaxlen
        self.ampminlen = ampminlen
        self.tmdiff = float(tmdiff)
        self.regionstart = regionstart
        self.regionend = regionend
        self.polyT = 8
    def print(self):
        s = ""
        s += "\n"
        s += "*" * 30 + "\n"
        s += "PrimerMaxLength:" + str(self.primermaxlen) + "\n"
        s += "PrimerMinLength:" + str(self.primerminlen) + "\n"
        s += "MaxTM:" + str(self.maxtm) + "\n"
        s += "MinTM:" + str(self.mintm) + "\n"
        s += "AmpMaxLen:" + str(self.ampmaxlen) + '\n'
        s += "AmpMinLen:" + str(self.ampminlen) + "\n"
        s += "TMDiff:" + str(self.tmdiff) + "\n"
        s += "RegionStart:" + str(self.regionstart) + "\n"
        s += "RegionEnd:" + str(self.regionend) + "\n"
        s += "*" * 30 + "\n"
        s += "\n"
        return s
    

