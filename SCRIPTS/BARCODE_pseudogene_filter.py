e#!/usr/bin/env python
# -*-coding: utf-8 -*

"""
Author : Benoit Penel
Date : 10/15/2024
version : 3.0
PYTHON3


run command : python3 BARCODE_pseudogene.py /path/to/DATA/RAW_ABUNDANCE_TABLE_fragment.csv 

This document filters the pseudogenous sequences of the input abundance tables

"""

import sys
import csv
import os 
import pandas as pd
import numpy 


Abundance_file=sys.argv[1]

nameprojet="AgriB07" #to modify put the name of your project
fragment="CR" # to modify CR or BB
outputdir="/home/lbenoit/Documents/Programmes/HAMI/DATA/RUN_AGRIB07/AgriB07_bar/"

DNA_length=int(325) #To modify according to CR (325) or BB (418)
reading_frame=int(1) #Check it before use it according to primer used
stop_codon=["TAA","TAG"] 


###OPEN FILES
csvfile=open(Abundance_file)
BARCODE = pd.read_csv(csvfile,delimiter="\t",header=0)



###############################################################
### FILTRE PSEUDOGENE (Length and stop codon) #################
###############################################################

seq=pd.Series(BARCODE.seed_sequence)
name=pd.Series(BARCODE.blast_taxonomy)


lenaccepted=[DNA_length-9,DNA_length-6,DNA_length-3,DNA_length,DNA_length+3,DNA_length+6,DNA_length+9] # Defines the range of DNA sequence lengths to be considered.  
accepted=[]
deleted=[]


print("Pseudogene research in ", fragment)

# CHECK INDEL IN EACH SEQUENCES
for i in range(len(seq)):
    if len(seq.iloc[i]) not in lenaccepted: #Keep position in the file of sequences who have undergone indels event(s) not a multiple of three (i.e. sequences with a different number of nucleotides than previously defined)
        deleted.append(i)
    else:

        #CHECK STOP CODON IN EACH SEQUENCES WITHOUT INDEL 
        catch = numpy.arange(reading_frame, len(seq.iloc[i]), 3) #Allows codons of three nucleotides to be read (indicated by the last argument = 3). The 1 can be changed to 0, 1 or 2 to shift the reading frame.
        stopCodonPositions = []
        bool=False
        for j in catch:
            codon = seq.iloc[i][j:j + 3]
            if codon in stop_codon:  # #Keep position in the file of sequences who have stop codon  
                deleted.append(i)
                bool=True
                break

        if bool==False: 
                accepted.append(i)          



pseudogene=BARCODE.iloc[deleted] #Keep pseudogenic OTU using position saved in deleted list 
realcol=BARCODE.iloc[accepted] #Keep non pseudogenic OTU using position saved in accepted list  

realcol.to_csv(os.path.join(outputdir+nameprojet+fragment+'_abundance_pseudogene_filtred.csv'),index_label=False,index=False,sep=';') #Abundance table without lines associated with pseudogenes


if not os.path.exists(os.path.join(outputdir+"pseudogenes/")):#Checks whether the "pseudogenes" directory exists in the file directory tree 
    os.makedirs(os.path.join(outputdir+"pseudogenes/"))

pseudogene.to_csv(os.path.join(outputdir+"pseudogenes/"+nameprojet+fragment+'_pseudogene_deleted.csv'),index_label=False,index=False,sep=';') #Abundance table with lines associated with pseudogenes

print("Pseudogene research in ", fragment, " over")

