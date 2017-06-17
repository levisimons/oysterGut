import os,re,sys
import gzip
from itertools import izip

leftRead = gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r')
rightRead = gzip.open('Undetermined_S0_L001_R2_001.fastq.gz','r')
rightIndex = gzip.open('Undetermined_S0_L001_I2_001.fastq.gz','r')
for leftRead, rightRead, rightIndex in izip(leftRead, rightRead, rightIndex):
    print "%s\t%s" % (leftRead.rstrip(), rightRead.rstrip(), rightIndex.rstrip())
