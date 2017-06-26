import os,re,sys
import gzip
from itertools import izip

lineNum=1

leftIndexList = ['TCAGC','GTATC','GCTAC','ACGCA',',CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT',
                 'GGCAG','TGATA','TGTGC']
#leftIndexList represents the forward indices
rightIndexList = ['CGTGA','ACATC','GCCTA','TGGTC','CACTG','ATTGG','GATCT','TCAAG','TGACA','GGACG','GCGGA','TTTCA','CCGGT','ATCGT','TGAGT','CGCCT']
#rightIndexList represents the reverse indices
leftFileSuffix = '_L001_R1_001.fastq'
rightFileSuffix = '_L001_R2_001.fastq'
leftIndex = 'Blank'
rightIndex = 'Blank'

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
#Note: the forward indices are on positions 4-9 (Counting from 0).  The 515F forward primer sequence 'GTGYCAGCMGCCGCGGTAA' is at positions 10-28.
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
#Note: the 926R primer sequence 'CCGYCAATTYMTTTRAGTTT' goes from positions 0-18.
indexRead = gzip.open('Undetermined_S0_L001_I1_001.fastq.gz','r')
#Note: The reverse indices are on positions 0-5.
for leftLine, rightLine, indexLine in izip(leftRead, rightRead, indexRead):
  if lineNum%4 == 1:
    leftIdentifier=str(leftLine)
    rightIdentifier=str(rightLine)
    #Read in the sequence identifier.  This the first of four lines.
  if lineNum%4 == 2:
    leftIndex=leftLine[4:9]
    rightIndex=indexLine[0:5]
    #Read in the forward and reverse indices.
    leftSequence=leftLine[28:]
    rightSequence=rightLine[18:]
    #Read in the forward and reverse sequences.
  if lineNum%4 == 3:
    leftQI = str(leftLine)
    rightQI = str(rightLine)
    #Somewhat redundant, but the third line of a fastq file is just '+'.
  if lineNum%4 == 0:
    leftQuality = leftLine[28:]
    rightQuality = rightLine[18:]
    if leftIndex in leftIndexList:
      if rightIndex in rightIndexList:
        #If the forward and reverse reads both have indices in the list then begin sorting the reads from the source fastq file into smaller fastq files by forward-reverse indices.
        leftIndexNum = leftIndexList.index(leftIndex)+1
        rightIndexNum = rightIndexList.index(rightIndex)+1
        leftFileName = 'F',str(leftIndexNum),'R',str(rightIndexNum),str(leftFileSuffix)
        leftFileName = ''.join(leftFileName)
        rightFileName = 'F',str(leftIndexNum),'R',str(rightIndexNum),str(rightFileSuffix)
        rightFileName = ''.join(rightFileName)
        leftOutput = open(leftFileName, 'a')
        leftOutputLine = leftIdentifier,leftSequence,leftQI,leftQuality
        leftOutputLine = ''.join(leftOutputLine)
        leftOutput.write(leftOutputLine)
        rightOutput = open(rightFileName, 'a')
        rightOutputLine = rightIdentifier,rightSequence,rightQI,rightQuality
        rightOutputLine = ''.join(rightOutputLine)
        rightOutput.write(rightOutputLine)
        leftIndex = 'Blank'
        rightIndex = 'Blank'
  lineNum=lineNum+1
