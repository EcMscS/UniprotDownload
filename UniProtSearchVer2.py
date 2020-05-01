import os
from UniprotObjectFasta import UniProt
import numpy as np
import requests

#Data will be in the following format: >organism | unitproID | description

#Test Case
sample1 = UniProt('sp', 'Q0P9C4', 'PGLK_CAMJE Protein glycosylation K OS=Campylobacter jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168) OX=192222 GN=pglK PE=1 SV=1', '')

#Reading data from fasta file format
path = "/Users/jeffreylai/Developer/Python/UniprotDownload/"
fileName = "Uniprot_Data_SampleFasta.txt"

try:
    file = open(path + fileName, 'r')
except IOError:
    print("An error was found. Either path is incorrect or file doesn't exist!")

#Extract line with > at the beginning of the line from file
extractedLineFromSearch = []
extractedProteinSequence = []

tempSeq = ""
seqNum = 1

for line in file:
    if line.startswith('>'):
        extractedLineFromSearch.append(line)
        if seqNum != 1:
            extractedProteinSequence.append(tempSeq)
            tempSeq = ""
        seqNum = seqNum + 1
    else:
        tempSeq = tempSeq + line.rstrip('\n')

extractedProteinSequence.append(tempSeq) #Get the last sequence
file.close

print("Number of Protein Sequences " + str(len(extractedProteinSequence)))

#number = 1
#for proteinSeq in extractedProteinSequence:
#    print(str(number) + " " + proteinSeq)
#    number = number + 1

#Create new UniProt Objects and make into array
uniprot_objects = []

for eachItem in extractedLineFromSearch:
    organism, uniprotID, description = eachItem.split('|')
    fullSequence = ''
    uniprotItem = UniProt(organism.strip('>'), uniprotID, description, fullSequence)
    uniprot_objects.append(uniprotItem)

index = 0
for eachSeq in extractedProteinSequence:
    uniprot_objects[index].fullSequence = eachSeq
    index = index + 1

for each in uniprot_objects:
    info = each.organism + " " + each.uniprotID + " " + each.description + " " + each.fullSequence

#sequenceListAsFasta = []
#or singleHttpRequest in uniprot_objects:
#    URL = UniProt.urlBase + singleHttpRequest.uniprotID + UniProt.urlFileType
#    r = requests.get(url=URL)
#    sequenceListAsFasta.append(r.text)

try:
    file = open(path + fileName, 'r')
except IOError:
    print("An error was found. Either path is incorrect or file doesn't exist!")

#Create list with only the second line of each fasta sequence and save as full sequence
extractedLine2FromFasta = []

index = 0
for eachOrganism in extractedLineFromSearch:
    searchWord = '\n'
    search = eachOrganism.find(searchWord)
    positionToExtract = search + 1
    sequenceEnd = len(eachOrganism) - 1
    fullSequence = eachOrganism[positionToExtract:positionToExtract + sequenceEnd].replace('\n','')
    extractedLine2FromFasta.append(fullSequence)
    print(str(index) + " " + uniprot_objects[index].fullSequence)
    uniprot_objects[index].fullSequence = fullSequence
    print(fullSequence)
    index = index + 1

file.close

print('Process Completed ; Created ' + str(len(uniprot_objects)) + ' uniprotObjects')
