import os,re,sys
import gzip
from itertools import izip

lineNum=1

leftIndexList = ['TCAGC','GTATC','GCTAC','ACGCA',',CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT',
                 'GGCAG','TGATA','TGTGC']

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
rightIndex = gzip.open('Undetermined_S0_L001_I1_001.fastq.gz','r')
for leftLine, rightLine, indexLine in izip(leftRead, rightRead, rightIndex):
    if lineNum%4 == 2:
        leftIndex=leftLine[4:9]
        if leftIndex in leftIndexList: primerNum = leftIndexList.index(leftIndex)
        print (leftIndex.rstrip(), primerNum, leftLine.rstrip(), rightLine.rstrip(), indexLine.rstrip())
    lineNum=lineNum+1
