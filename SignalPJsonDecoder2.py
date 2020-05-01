import json
import os
from UniprotObjectFasta import UniProt
import pandas as pd
import requests

#Code to parse JSON
def extract_values(obj, key):
    """Pull all values of specified key from nested JSON."""
    arr = []

    def extract(obj, arr, key):
        """Recursively search for values of key in JSON tree."""
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, (dict, list)):
                    extract(v, arr, key)
                elif k == key:
                    arr.append(v)
        elif isinstance(obj, list):
            for item in obj:
                extract(item, arr, key)
        return arr

    results = extract(obj, arr, key)
    return results

#Change Variables here before running code
path = "/Users/jeffreylai/Developer/Python/UniprotDownload/"
fileName = "output_1623_Seq.json" #SignalP JSON Result Summary

fileNameOriginal = "Uniprot_Data_for_SignalP.txt" #Uniprot Fasta file submitted to SignalP

#Getting Fasta Sequences
#--------------
#Read from File
#Identifiers: Organism, UniprotID, Description, Full Sequence
try:
    file = open(path + fileNameOriginal, 'r')
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
#for singleHttpRequest in uniprot_objects:
#     URL = UniProt.urlBase + singleHttpRequest.uniprotID + UniProt.urlFileType
#     r = requests.get(url=URL)
#     print(r.text)
#     sequenceListAsFasta.append(r.text)

#Create list with only the second line of each fasta sequence and save as full sequence
#extractedLine2FromFasta = []

#sequenceListAsFasta -> fileNameOriginal
#index = 0
#or eachOrganism in fileNameOriginal:
#    searchWord = '\n'
#    search = eachOrganism.find(searchWord)
#    positionToExtract = search + 1
#    sequenceEnd = len(eachOrganism) - 1
#    fullSequence = eachOrganism[positionToExtract:positionToExtract + sequenceEnd].replace('\n','')
#    extractedLine2FromFasta.append(fullSequence)
#    print("Full Seq: " + fullSequence)
#    uniprot_objects[index].fullSequence = fullSequence
#    index = index + 1

#---------------------------------------

with open(path + fileName, "r") as read_file:
    data = json.load(read_file)

#Prints all sequences from JSON from top level
print(data['SEQUENCES'])

#Parse JSON Function
organism = extract_values(data, 'Name')
cleavageSite = extract_values(data, 'CS_pos')
predictionSignalPeptide = extract_values(data, 'Prediction')
print(organism)
print(cleavageSite)
print(predictionSignalPeptide)

#Create 2 DataFrames
organismList = {'organism': organism}
cleavageSiteList = {'cleavage_site': cleavageSite}

organismDF = pd.DataFrame(organismList)
cleavageSiteDF = pd.DataFrame(cleavageSiteList)
predictionListDF = pd.DataFrame(predictionSignalPeptide)

combinedDF1 = pd.DataFrame.join(organismDF, cleavageSiteDF)
combinedFinalDF = pd.DataFrame.join(combinedDF1, predictionListDF)
print(combinedFinalDF)

#column_names = ['Organism', 'cleavage_site']
#column_datatype = ['string', 'string']

#schema_dict = dict(zip(column_names, column_datatype))
#print(schema_dict)

#Find pos.

#Find first case where there is a cleavage site
correctIndex = 0
for eachItem in cleavageSite:
    if len(eachItem) != 0:
        correctIndex = correctIndex
        break
    correctIndex = correctIndex + 1


searchWord = 'pos.'
search = cleavageSite[correctIndex].find(searchWord)
positionToExtract = search + 5
sampleCleavageSite = cleavageSite[0][positionToExtract: positionToExtract+2]
#print("SAMPLE: " + sampleCleavageSite)

#Extract all positions from cleavageSiteList
#With position number obtained from above
extractedCleavageSitePositionList = []

for eachItem in cleavageSite:
    if len(eachItem) != 0:
        #print("Cleavage: " + eachItem)
        sampleCleavageSite = eachItem[positionToExtract: positionToExtract+2]
        #print("Site: " + sampleCleavageSite)
        extractedCleavageSitePositionList.append(sampleCleavageSite)
    else:
        extractedCleavageSitePositionList.append(0)

#Create DateFrame of Cleavage cleavage sites
cleavagePositionNumber = {'Position':extractedCleavageSitePositionList}
cleavagePositionNumberDF = pd.DataFrame(cleavagePositionNumber)
combinedWithPositionNumber = pd.DataFrame.join(combinedDF1, cleavagePositionNumberDF)
print(combinedWithPositionNumber)

#Add fullSequence to dataFrame
fullSequenceOfOrganism = {'Full_Sequence':extractedProteinSequence}
fullSequenceOfOrganismDF = pd.DataFrame(fullSequenceOfOrganism)
addedFullSequence = pd.DataFrame.join(combinedWithPositionNumber, fullSequenceOfOrganismDF)
print(addedFullSequence)

#Cut sequences according to cleavage site
#Use extractedLine2FromFasta to create new signal peptide sequence
signalPeptideSequence = []

#PROBLEM HERE
# 3 is position and 4 is Full Sequence
for eachOrganism in addedFullSequence.itertuples():
    site = int(eachOrganism[3])
    signalPeptideSequence.append(eachOrganism[4][0:site])

#Add signal peptide sequence to dataFrame
#signalPeptideSequenceOfOrganism = {'Signal_Sequence':signalPeptideSequence}
#signalPeptideSequenceDF = pd.DataFrame(signalPeptideSequenceOfOrganism)
#addedCutSequence = pd.DataFrame.join(addedFullSequence, signalPeptideSequenceDF)

#Add signal peptide sequence to dataFrame
signalPeptideSequenceOfOrganism = {'Signal_Sequence':signalPeptideSequence}
signalPeptideSequenceDF = pd.DataFrame(signalPeptideSequenceOfOrganism)
addedSignalPeptideSequence = pd.DataFrame.join(addedFullSequence, signalPeptideSequenceDF)

#Determine if protein is Lipoprotein
#Check if Predition has the lipoprotein in the results
isLipoproteinList = []

for eachOrganism in predictionSignalPeptide:
    findLipoprotein = eachOrganism.find('Lipoprotein')
    #print("Prediction Peptide: " + eachOrganism)
    if findLipoprotein != -1:
        isLipoproteinList.append('Yes')
    else:
        isLipoproteinList.append('No')

#Create DateFrame on whether Organism is a lipoprotein
isLipoproteinResults = {'Is Lipoprotein': isLipoproteinList}
isLipoproteinResultsDF = pd.DataFrame(isLipoproteinResults)
combinedWithIsLipoprotein = pd.DataFrame.join(addedSignalPeptideSequence, isLipoproteinResultsDF)
print(combinedWithIsLipoprotein)

#Export to Excel
combinedWithIsLipoprotein.to_excel(r'/Users/jeffreylai/Developer/Python/UniprotDownload/signalSequenceResults1623.xlsx')

#create new dataframe and add join together

print("Process Completed ; SignalP JSON Results Parsed")
