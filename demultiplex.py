import os,re,sys
import gzip

with gzip.open('Undetermined_S0_L001_R1_001.fastq.gz','r') as rightIndex:        
            for line in rightIndex:        
                print('got line', line)
