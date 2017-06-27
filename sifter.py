import os,re,sys
from subprocess import call

inputFile = open('unused.txt','rU')

for line in inputFile:
  line=line.rstrip('\r\n')
  command1 = 'mv ',line,'_L001_R1_001.fastq unused'
  command1 = ''.join(command1)
  command2 = 'mv ',line,'_L001_R2_001.fastq unused'
  command2 = ''.join(command2)
  print command1
  print command2
  call(command1)
  call(command2)
