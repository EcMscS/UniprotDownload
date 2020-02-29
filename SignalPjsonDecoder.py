import json
import os
from UniProtObject import UniProt
import pandas as pd
import requests

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

path = "/Users/jeffreylai/Developer/Python/Uniprot/"
fileName = "output.json"

#Getting Fasta Sequences
#--------------
#Read from File
#Identifiers: Organism, Locus Tag, UniProtID, isLipoprotein
fileNameOriginal = "SampleData.txt"
file = open(path + fileNameOriginal, 'rt')

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

sequenceListAsFasta = []
for singleHttpRequest in uniprot_objects:
    URL = UniProt.urlBase + singleHttpRequest.uniprotID + UniProt.urlFileType
    r = requests.get(url=URL)
    sequenceListAsFasta.append(r.text)

#Create list with only the second line of each fasta sequence
extractedLine2FromFasta = []

for eachOrganism in sequenceListAsFasta:
    searchWord = '\n'
    search = eachOrganism.find(searchWord)
    positionToExtract = search + 1
    sequenceEnd = len(eachOrganism) - 1
    fullSequence = eachOrganism[positionToExtract:positionToExtract + sequenceEnd].replace('\n','')
    extractedLine2FromFasta.append(fullSequence)

print(extractedLine2FromFasta)
#--------------


with open(path + fileName, "r") as read_file:
    data = json.load(read_file)

#Prints all sequences from JSON from top level
print(data['SEQUENCES'])

organism = extract_values(data, 'Name')
cleavageSite = extract_values(data, 'CS_pos')
print(organism)
print(cleavageSite)

#Create 2 DataFrames
organismList = {'organism': organism}
cleavageSiteList = {'cleavage_site': cleavageSite}

organismDF = pd.DataFrame(organismList)
cleavageSiteDF = pd.DataFrame(cleavageSiteList)

combined = pd.DataFrame.join(organismDF, cleavageSiteDF)
print(combined)

column_names = ['Organism', 'cleavage_site']
column_datatype = ['string', 'string']

#schema_dict = dict(zip(column_names, column_datatype))
#print(schema_dict)

#Find pos.
searchWord = 'pos.'
search = cleavageSite[0].find(searchWord)
positionToExtract = search + 5
sampleCleavageSite = cleavageSite[0][positionToExtract: positionToExtract+2]
print(sampleCleavageSite)

#Extract all positions from cleavageSiteList
#With position number obtained from above
extractedCleavageSitePositionList = []

for eachItem in cleavageSite:
    print(eachItem)
    sampleCleavageSite = eachItem[positionToExtract: positionToExtract+2]
    extractedCleavageSitePositionList.append(sampleCleavageSite)

#Create DateFrame of Cleavage cleavage sites
cleavagePositionNumber = {'Position':extractedCleavageSitePositionList}
cleavagePositionNumberDF = pd.DataFrame(cleavagePositionNumber)
combinedWithPositionNumber = pd.DataFrame.join(combined, cleavagePositionNumberDF)
print(combinedWithPositionNumber)

#Add fullSequence to dataFrame
fullSequenceOfOrganism = {'Full_Sequence':extractedLine2FromFasta}
fullSequenceOfOrganismDF = pd.DataFrame(fullSequenceOfOrganism)
addedFullSequence = pd.DataFrame.join(combinedWithPositionNumber, fullSequenceOfOrganismDF)
#print(addedFullSequence)

#Cut sequences according to cleavage site
#Use extractedLine2FromFasta to create new signal peptide sequence
signalPeptideSequence = []

# 3 is position and 4 is Full Sequence
for eachOrganism in addedFullSequence.itertuples():
    site = int(eachOrganism[3])
    signalPeptideSequence.append(eachOrganism[4][0:site])

#Add signal peptide sequence to dataFrame
signalPeptideSequenceOfOrganism = {'Signal_Sequence':signalPeptideSequence}
signalPeptideSequenceDF = pd.DataFrame(signalPeptideSequenceOfOrganism)
finalDF = pd.DataFrame.join(addedFullSequence, signalPeptideSequenceDF)

print(finalDF)

#Export to Excel
finalDF.to_excel(r'/Users/jeffreylai/Developer/Python/Uniprot/signalSequence.xlsx')



#create new dataframe and add join together
