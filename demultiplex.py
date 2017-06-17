import os,re,sys
import gzip
from itertools import izip

lineNum=1

leftIndexList = ['TCAGC','GTATC','GCTAC','ACGCA',',CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT',
                 'GGCAG','TGATA','TGTGC']
rightIndexList = ['CGTGA','ACATC','GCCTA','TGGTC','CACTG','ATTGG','GATCT','TCAAG','TGACA','GGACG','GCGGA','TTTCA','CCGGT','ATCGT','TGAGT','CGCCT']
leftFileSuffix = '_L001_R1_001.fastq'
rightFileSuffix = '_L001_R2_001.fastq'
match = False

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
indexRead = gzip.open('Undetermined_S0_L001_I1_001.fastq.gz','r')
for leftLine, rightLine, indexLine in izip(leftRead, rightRead, indexRead):
    if lineNum%4 == 1:
      leftIdentifier=str(leftLine)
      rightIdentifier=str(rightLine)
    if lineNum%4 == 2:
        leftIndex=leftLine[4:9]
        rightIndex=indexLine[0:5]
        if leftIndex in leftIndexList:
          if rightIndex in rightIndexList:
            match = True
            leftSequence=leftLine[9:]
            rightSequence=rightLine[9:]
            leftIndexNum = leftIndexList.index(leftIndex)+1
            rightIndexNum = rightIndexList.index(rightIndex)+1
            #print 'Left ', leftIndexNum, leftIndex, 'Right ', rightIndexNum, rightIndex
            leftFileName = 'F',str(leftIndexNum),'_R',str(rightIndexNum),str(leftFileSuffix)
            leftFileName = ''.join(leftFileName)
            rightFileName = 'F',str(leftIndexNum),'_R',str(rightIndexNum),str(rightFileSuffix)
            rightFileName = ''.join(rightFileName)
            #print leftFileName, rightFileName
          else: match = False
        else: match = False
    if lineNum%4 == 3:
      leftQI = str(leftLine)
      rightQI = str(rightLine)
    if lineNum%4 == 0:
      leftQuality = str(leftLine)
      rightQuality = str(rightLine)
    if match == True:
      leftOutput = open(leftFileName, 'a')
      leftOutputLine = leftIdentifier,leftSequence,leftQI,leftQuality
      leftOutputLine = ''.join(leftOutputLine)
      leftOutput.write(leftOutputLine)
      rightOutput = open(rightFileName, 'a')
      rightOutputLine = rightIdentifier,rightSequence,rightQI,rightQuality
      rightOutputLine = ''.join(rightOutputLine)
      rightOutput.write(rightOutputLine)
    lineNum=lineNum+1
