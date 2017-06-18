import os,re,sys
import gzip
from itertools import izip

lineNum=1

leftIndexList = ['TCAGC','GTATC','GCTAC','ACGCA',',CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT',
                 'GGCAG','TGATA','TGTGC']
rightIndexList = ['CGTGA','ACATC','GCCTA','TGGTC','CACTG','ATTGG','GATCT','TCAAG','TGACA','GGACG','GCGGA','TTTCA','CCGGT','ATCGT','TGAGT','CGCCT']
leftFileSuffix = '_L001_R1_001.fastq'
rightFileSuffix = '_L001_R2_001.fastq'
leftIndex = 'Blank'
rightIndex = 'Blank'

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
indexRead = gzip.open('Undetermined_S0_L001_I1_001.fastq.gz','r')
for leftLine, rightLine, indexLine in izip(leftRead, rightRead, indexRead):
  leftIndex = 'Blank'
  rightIndex = 'Blank'
  if lineNum%4 == 1:
    leftIdentifier=str(leftLine)
    rightIdentifier=str(rightLine)
  if lineNum%4 == 2:
    leftIndex=leftLine[4:9]
    rightIndex=indexLine[0:5]
  if lineNum%4 == 3:
    leftQI = str(leftLine)
    rightQI = str(rightLine)
  if lineNum%4 == 0:
    leftQuality = str(leftLine)
    rightQuality = str(rightLine)
    if leftIndex in leftIndexList:
      if rightIndex in rightIndexList:
        leftSequence=leftLine[9:]
        rightSequence=rightLine[0:]
        leftIndexNum = leftIndexList.index(leftIndex)+1
        rightIndexNum = rightIndexList.index(rightIndex)+1
        #print 'Left ', leftIndexNum, leftIndex, 'Right ', rightIndexNum, rightIndex
        leftFileName = 'F',str(leftIndexNum),'_R',str(rightIndexNum),str(leftFileSuffix)
        leftFileName = ''.join(leftFileName)
        rightFileName = 'F',str(leftIndexNum),'_R',str(rightIndexNum),str(rightFileSuffix)
        rightFileName = ''.join(rightFileName)
        #print leftFileName, rightFileName
        #leftOutput = open(leftFileName, 'a')
        leftOutputLine = leftIdentifier,leftSequence,leftQI,leftQuality
        leftOutputLine = ''.join(leftOutputLine)
        #leftOutput.write(leftOutputLine)
        print leftFileName,leftOutputLine
        rightOutput = open(rightFileName, 'a')
        rightOutputLine = rightIdentifier,rightSequence,rightQI,rightQuality
        rightOutputLine = ''.join(rightOutputLine)
        #rightOutput.write(rightOutputLine)
        print rightFileName,rightOutputLine
  lineNum=lineNum+1
