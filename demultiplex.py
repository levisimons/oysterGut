import os,re,sys
import gzip
from itertools import izip

lineNum=1

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
rightIndex = gzip.open('Undetermined_S0_L001_I2_001.fastq.gz','r')
for leftLine, rightLine, indexLine in izip(leftRead, rightRead, rightIndex):
    if lineNum%4 == 2:
        leftIndex=leftLine[4:9]
        print (leftIndex.rstrip(), leftLine.rstrip(), rightLine.rstrip(), indexLine.rstrip())
    lineNum=lineNum+1
