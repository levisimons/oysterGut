import os,re,sys
import gzip
from itertools import izip

lineNum=1

leftIndexList = ['TCAGC','GTATC','GCTAC','ACGCA',',CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT',
                 'GGCAG','TGATA','TGTGC']
rightIndexList = ['CGTGA','ACATC','GCCTA','TGGTC','CACTG','ATTGG','GATCT','TCAAG','TGACA','GGACG','GCGGA','TTTCA','CCGGT','ATCGT','TGAGT','CGCCT']
leftFileSuffix = 'L001_R1_001.fastq'
rightFileSuffix = 'L001_R2_001.fastq'

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
indexRead = gzip.open('Undetermined_S0_L001_I1_001.fastq.gz','r')
for leftLine, rightLine, indexLine in izip(leftRead, rightRead, indexRead):
    if lineNum%4 == 2:
        leftIndex=leftLine[4:9]
        rightIndex=indexLine[0:5]
        if leftIndex in leftIndexList:
          if rightIndex in rightIndexList:
            leftIndexNum = leftIndexList.index(leftIndex)+1
            rightIndexNum = rightIndexList.index(rightIndex)+1
            print 'Left ', leftIndexNum, leftIndex, 'Right ', rightIndexNum, rightIndex
            leftFileName = ''.join('F',leftIndexNum,'R',rightIndexNum,leftFileSuffix)
            rightFileName = ''.join('F',leftIndexNum,'R',rightIndexNum,rightFileSuffix)
            print leftFileName, rightFileName
    lineNum=lineNum+1
