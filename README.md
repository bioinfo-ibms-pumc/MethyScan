# MethyScan: A tool for methylation specific PCR primer design and evaluation

MethyScan can simultaneously complete the design of two MSP primers and the nested primer for increased sensitivity. In addition, by converting the genome sequence and integrating the next-generation sequencing tool [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml), MethyScan can evaluate the non-specific amplification information of the primers on the genome sequence. 

MethyScan is maintained by Yinghao Cao [yhcao@ibms.pumc.edu.cn].

## Download and Installation
```
git clone https://github.com/bioinfo-ibms-pumc/MethyScan.git
```

## Dependency package installation:
```
pip3 install biopython matplotlib primer3-py
```
  MethyScan also needs other NGS softwares: [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml), [BEDTools](https://bedtools.readthedocs.io/en/latest/), [Samtools](http://samtools.sourceforge.net/)


## Command Lines

```  
usage: MethyScan.py

Program: MethyScan
Version: 1.0
Contact: Yinghao Cao <yhcao@ibms.pumc.edu.cn>
	

positional arguments:
  {designprimer,formatdb,searchdb}
                        command
    designprimer        Design primer pairs for target sequence.
    formatdb            Create methylated databse.
    searchdb            Search primer pairs against methylated database.

optional arguments:
  -h, --help            show this help message and exit
  
usage A: MethyScan.py designprimer [-h] -i INPUTFILE [-o OUTPUTFILE] 
                        [--minlength MINLENGTH] [--maxlength MAXLENGTH]
                        [--mintm MINTM] [--maxtm MAXTM] [-m MAXNUM]
                        [--tmdiff TMDIFF]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputfile INPUTFILE
                        Target sequence file in fasta format.
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        Output file.
  --minlength MINLENGTH
                        Primer minimum length. [Default:20]
  --maxlength MAXLENGTH
                        Primer maximum length. [Default:25]
  --mintm MINTM         Primer minimum TM. [Default:50]
  --maxtm MAXTM         Primer maximum TM. [Default:60]
  -m MAXNUM, --maxnum MAXNUM
                        Maximum primer pairs for one CpG island. [Default:3]
  --tmdiff TMDIFF       Maximum TM difference in primer pairs. [Default:5]

usage B: MethyScan.py formatdb [-h] -g GENOMEFILE [-o OUTPUTFILE] [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMEFILE, --genomefile GENOMEFILE
                        Database file in fasta format.
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        Output database name.
  -t THREADS, --threads THREADS
                        Number of threads used for building index in bowtie. [Default:4]

usage C: MethyScan.py searchdb [-h] -i INPUTFILE -d DATABASE -o OUTPUTFILE 
                        [--minlength MINLENGTH] [--maxlength MAXLENGTH]
                        [-c MISMATCH]  [-k SEARCHTYPE] [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputfile INPUTFILE
                        Input primer sequence file in tab format likes "amp1 p5seq p3seq")
  -d DATABASE, --database DATABASE
                        Methylated database filename prefix.
  -c MISMATCH, --mismatch MISMATCH
                        Maximum mismatch in one primer sequence, (<=3). [Default:3]
  --minlength MINLENGTH
                        Minimum amplicon length considered in mapping. [default:100]
  --maxlength MAXLENGTH
                        Maximum amplicon length considered in mapping. [default:500]
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        Output file
  -k SEARCHTYPE, --searchtype SEARCHTYPE
                        Search against methylated database with specific strand. Type: ALL,[MCT],[MGA],[UCT],[UGA]. [Default:ALL]
  -t THREADS, --threads THREADS
                        Number of threads used for searching. [default: 1]

```
## Examples
```
1. Design primers for target sequences.
   python3 MethyScan.py designprimer -i tfpi2.fa -o mytest
   
2. Genome sequence conversion and indexing.
   python3 MethyScan.py formatdb -i Homo_sapiens_assembly38.fasta -o mydb

3. Search primer pairs against methylated genome sequences.
   python3 MethyScan.py searchdb -i mytest.txt -d mydb -t 10 -o primer_res.txt -k ALL
  
```
## Citing

If you use MethyScan for your research, please kindly cite the following paper:

曹英豪.MethyScan：一种甲基化特异性PCR引物设计及评估工具[J].生物化学与生物物理进展,2021,48(06):677-687.
CAO Ying-Hao.MethyScan：A Tool for Methylation Specific PCR Primer Design and Evaluation[J].Progress in Biochemistry and Biophysics,2021,48(06):677-687. 
DOI:[https://doi.org/10.16476/j.pibb.2020.0413](https://doi.org/10.3389/fgene.2020.00490)
