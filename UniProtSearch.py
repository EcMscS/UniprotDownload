import os
from UniProtObject import UniProt
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import requests

#Test UniProt
test1 = UniProt('N_meningitidis', 'NMB1946', 'Q7DD63' ,'yes')

#Read from File
#Identifiers: Organism, Locus Tag, UniProtID, isLipoprotein
path = "/Users/jeffreylai/Developer/Python/Uniprot/"
fileName = "SampleData.txt"
file = open(path + fileName, 'rt')

#Extract information from file
listOfEachSearch = []
for line in file:
    listOfEachSearch.append(line)
file.close

#Clean up strip line of extraneous info
separated = []
for line in listOfEachSearch:
    separated.append(line.strip())

#Create new UniProt Objects and make into array
uniprot_objects = []

for eachItem in separated:
    organism, locus_tag, id, lipoprotein = eachItem.split(',')
    uniprot = UniProt(organism.strip('\''), locus_tag.strip('\''), id.strip('\''), lipoprotein.strip('\''))
    uniprot_objects.append(uniprot)

for each in uniprot_objects:
    info = each.organism + " " + each.locusTag + " " + each.uniprotID + " " + each.isLipoprotein
    print(info)

sequenceListAsFasta = []

#Retrieve fasta seqence from UniProt
for singleHttpRequest in uniprot_objects:
    URL = UniProt.urlBase + singleHttpRequest.uniprotID + UniProt.urlFileType
    r = requests.get(url=URL)
    sequenceListAsFasta.append(r.text)

saveFile = open(path + "dnaResult.txt", "w")

for each in sequenceListAsFasta:
    saveFile.write(each)

saveFile.close

print("Process Completed")
