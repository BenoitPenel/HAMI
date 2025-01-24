#!/usr/bin/env python
# -*-coding: utf-8 -*

'''

Author : Benoit Penel
Date : 28/08/2023
version : 1.0
PYTHON3

The script is to discriminate Metabarcoding, Barcoding and control samples which were sequenced at the same time in a context of HAMI FRAMEWORK.
In practical, it separates  Metabarcoding, Barcoding and control datas from the input file which were - previously produce with frogs- by using the sample name nomenclature .
Your samples thus need to be discrimated on the basis of their name. Please use alphabectic prefix and number for it e.g : CMEY0001 
also pipeline is developped to handle duplicate of samples. Discrimination between duplicate will be done using suffix : e.g : CMEY0001A / CMEY0001B
    
According to the input file, the script will produce 3 output files : two of them are subset  of input file (Metabarcoding samples and Barcodign samples). The latter is METADATA file associated to metabarcoding datas
    
Note : If no barcoding sample, the Barcoding subfile will be empty.

'''

###### PACKAGES :

import sys
import pandas as pd
import os 
import re
import csv

###### PARAM YOUR PROJECT :
#Directory path
outputdir= sys.argv[2]
#Name file
name=sys.argv[3]
fragment=sys.argv[4]
#Parameter necessary to code different type of samples
prefixMETA = sys.argv[5]
prefixBARC = sys.argv[6]
prefixCONTROL = sys.argv[7]
duplicat =sys.argv[8]
duplicatformule=duplicat.replace("/","|")



###############
#### BEGIN ####
###############


###OPEN FILES
csvfile = open(sys.argv[1])
Occurencefile = pd.read_csv(csvfile, delimiter="\t", header=0)

colname = list(Occurencefile.columns.values)#Create an object using column headers 

first=colname[0:11] #Keep the 11th first colum of the dataframe previously produce by FROGS. It contains information about taxonomic affiliations, identity percentages, sequences etc...



### FILTERING SAMPLES ACCORDING TO THEIR ID_NAME, WHICH INDICATES WHETHER THEY ARE BARCODING, METABARCODING, CONTROL SAMPLES 
BARCODE = []
METABARCODE = []
CONTROL = []

for i in colname[11:]:   # 11th first colum are already keep in "first" variable
    if re.search(prefixBARC, i): #search barcoding samples using prefix nomenclature 
        BARCODE.append(i)
    else : 
        if re.search(prefixMETA, i): #search metabarcoding samples using prefix nomenclature
            METABARCODE.append(i)
        if re.search(prefixCONTROL+r'(I|E|P)', i, re.IGNORECASE): #search control samples using prefix nomenclature
            CONTROL.append(i)

### PRODUCE SUBSET OF THE INPUT FILE  :  

## METABARCODING
#Metabarcoding subfile : 
METABARCODE=first+METABARCODE+CONTROL
METAfile= Occurencefile[METABARCODE]
# Do the read sum to update "observation_sum" variable of the subfile :
sum_by_row = METAfile.iloc[:, 11:].sum(axis=1) 
METAfile.loc[:,"observation_sum"]=sum_by_row
METAfile= METAfile[METAfile["observation_sum"] != 0]
#Create metabarcoding file : 
if not os.path.exists(os.path.join(outputdir,"METABARCODING")):#Checks whether the "Metabarcoding" directory exust in the file directory tree 
    os.makedirs(os.path.join(outputdir,"METABARCODING"))
METAfile.to_csv(os.path.join(outputdir,"METABARCODING/"+name+fragment+"_abundance_raw_data_METABARCODING.tsv"), index=False,sep="\t")


##BARCODING
#barcoding subfile : 
if len(BARCODE) != 0: #If barcoding data 

    BARCODE=first+BARCODE+CONTROL
    BARCOfile= Occurencefile[BARCODE]
    # Do the read sum to update "observation_sum" variable of the subfile :
    sum_by_row = BARCOfile.iloc[:, 11:].sum(axis=1)
    BARCOfile.loc[:,"observation_sum"]=sum_by_row
    BARCOfile = BARCOfile[BARCOfile["observation_sum"]!= 0]
    #Create barcoding file : 

    if not os.path.exists(os.path.join(outputdir,"BARCODE")): #Checks whether the "barcoding" directory exust in the file directory tree 
        os.makedirs(os.path.join(outputdir,"BARCODE"))
    BARCOfile.to_csv(os.path.join(outputdir,"BARCODE/"+name+fragment+"_abundance_raw_data_BARCODING.tsv"), index=False,sep="\t")

else: # If no barcoding data 
    if not os.path.exists(os.path.join(outputdir,"BARCODE")):  #Checks whether the "barcoding" directory exust in the file directory tree 
        os.makedirs(os.path.join(outputdir,"BARCODE"))
    Null=os.path.join(outputdir,"BARCODE/"+name+fragment+"_abundance_raw_data_BARCODING.tsv")
    with open(Null, mode='w', newline='') as fichier_tsv: #Create an empty barcoding subfile
        writer = csv.writer(fichier_tsv)
        writer.writerow("NO BARCODING")


### PRODUCE METADATA for metabarcoding subfile : 

with open(os.path.join(outputdir,"METABARCODING/"+name+fragment+"_Metadata_METABARCODING.csv"),'w',newline='') as csvfile:
    csvwriter=csv.writer(csvfile,delimiter=";")
    header = ["observation_name", "rep", "control", "biological_unit", "project"]
    csvwriter.writerow(header)
    print(METABARCODE)
    for i in METABARCODE: #DUPLICATE ONLY (TRIPLICATE NEED CHANGE)
        print(prefixMETA,i)
        if re.search(prefixMETA, i) and  re.search(r'[-_]{}$'.format(duplicat[0]), i):
            row = [i, "1", "no", i[:-2], name]
            csvwriter.writerow(row)
        elif re.search(prefixMETA,i) and  re.search(r'[-_]{}$'.format(duplicat[-1]), i):
            row = [i, "2", "no", i[:-2], name]
            csvwriter.writerow(row)
        elif re.search(r'^NC(I|E|P)', i, re.IGNORECASE) and  re.search(r'[-_]{}$'.format(duplicat[0]), i):
            row= [i, "1", "negative", i[:-2], name]
            csvwriter.writerow(row)
        elif re.search(r'^NC(I|E|P)', i, re.IGNORECASE) and  re.search(r'[-_]{}$'.format(duplicat[-1]), i):
            row= [i, "2", "negative", i[:-2], name]
            csvwriter.writerow(row)
