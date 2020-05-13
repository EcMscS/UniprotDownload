#Modules needed for Python program to work
# coding=utf-8
import json
import os
import pandas as pd

from UniprotObjectFasta import UniProt #Check fro UniprtoObjectFasta.py file

#Function to parse JSON
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

#------------------------------------------------------------
####IMPORTANT####: These variables must be changed for the program to work
#Change Variables here before running code
path = "/Users/jeffreylai/Developer/Python/UniprotDownload/"

#These 2 files must be present in the above directory
fileName = "output_1623_Seq.json" #SignalP JSON Result Summary
fileNameOriginal = "Uniprot_Data_for_SignalP.txt" #Uniprot Fasta (Protein) file submitted to SignalP

exportedExcelFileName = "signalSequenceResults1623.xlsx" #Excel file name that contains all the extracted data

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

#Initialize variable for for loop
tempSeq = ""
seqNum = 1

#Extract information from fasta file, separate into: Organism, UniprotID, Description, and Full Protein Sequence
#extractedLineFromSearch array contains Organism | UniprotID | description
#extractedProteinSequence array contains each organisms full protein sequence
for line in file:
    if line.startswith('>'):
        extractedLineFromSearch.append(line) #Extract line with > with starting character
        if seqNum != 1:
            extractedProteinSequence.append(tempSeq) #Get the first sequence
            tempSeq = "" #Clear sequence for next iteration
        seqNum = seqNum + 1 #Increment current protein sequence
    else:
        tempSeq = tempSeq + line.rstrip('\n') #Append sequence with return carriage stripped

extractedProteinSequence.append(tempSeq) #Get the last sequence after the last for loop
file.close #Close file after operations completed

#Create new UniProt Objects and store as an array
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


#Separate into array for dataframes for later use
fastaOrganism = []
fastaUniprotID = []
fastaDescription = []
fastaFullSequence = []

for each in uniprot_objects:
    info = each.organism + " " + each.uniprotID + " " + each.description + " " + each.fullSequence
    fastaOrganism.append(each.organism)
    fastaUniprotID.append(each.uniprotID)
    fastaDescription.append(each.description)
    fastaFullSequence.append(each.fullSequence)

#---------------------------------------

with open(path + fileName, "r") as read_file:
    data = json.load(read_file)

#Testing JSON Parsing ; Prints all sequences from JSON from top level
#print(data['SEQUENCES'])

#Include fasta data
#Organism name, UniprotID and Description
fastaOrganismList = {'Organism_Fasta': fastaOrganism}
fastaUniprotIDList = {'UniprotID_Fasta': fastaUniprotID}
fastaDescriptionList = {'Description_Fasta': fastaDescription}
fastaFullSequenceList = {'Full_Sequence_Fasta': fastaFullSequence}

#Create Dataframes
organismFastaDF = pd.DataFrame(fastaOrganismList)
uniprotIDFastaDF = pd.DataFrame(fastaUniprotIDList)
descriptionFastaDF = pd.DataFrame(fastaDescriptionList)
fastaFullSequenceDF = pd.DataFrame(fastaFullSequenceList)

#Combine DataFrames
combinedFastaDF1 = pd.DataFrame.join(uniprotIDFastaDF, organismFastaDF)
combinedFastaDF2 = pd.DataFrame.join(combinedFastaDF1, descriptionFastaDF)
combinedFastaDF3 = pd.DataFrame.join(combinedFastaDF2, fastaFullSequenceDF)

fasta_Keys = combinedFastaDF3.keys()
keyFasta = fasta_Keys[0] #First element is UniProtID
combinedFastaDF3.sort_values(by=keyFasta, inplace=True)

#Sorted Full Sequence from FASTA File
fastaSortedSequence = combinedFastaDF3[['Full_Sequence_Fasta']]

#------------------------------------------------------------------------------
#Parse JSON Function from SignalP JSON Result
organism = extract_values(data, 'Name')
cleavageSite = extract_values(data, 'CS_pos')
predictionSignalPeptide = extract_values(data, 'Prediction')

#Extract UniprotID from Organism Name
#Sample Data: 'sp_P45491_Y984_CAMJE'
#After first _ is the UniProtID, ID is always 6 characters
#first, third, and fourth are not used.
organismUniProtID = []
for eachValue in organism:
    first, id, third, fourth = eachValue.split('_')
    organismUniProtID.append(id)

#Create DataFrame for UniProtID for SignalP Output
organismUniProtIDList = {'UniProtID': organismUniProtID}
organismUniProtIDDF = pd.DataFrame(organismUniProtIDList)

#Create 2 DataFrames
organismList = {'organism': organism}
cleavageSiteList = {'cleavage_site': cleavageSite}

organismDF = pd.DataFrame(organismList)
cleavageSiteDF = pd.DataFrame(cleavageSiteList)
predictionListDF = pd.DataFrame(predictionSignalPeptide)

combinedDF1 = pd.DataFrame.join(organismDF, cleavageSiteDF)
#combinedFinalDF = pd.DataFrame.join(combinedDF1, predictionListDF)

#Find first case where there is a cleavage site from sample
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

#Extract all positions from cleavageSiteList
#With position number obtained from above
extractedCleavageSitePositionList = []

for eachItem in cleavageSite:
    if len(eachItem) != 0:
        sampleCleavageSite = eachItem[positionToExtract: positionToExtract+2]
        extractedCleavageSitePositionList.append(sampleCleavageSite)
    else:
        extractedCleavageSitePositionList.append(0)

#Create DateFrame of Cleavage cleavage sites
cleavagePositionNumber = {'Position':extractedCleavageSitePositionList}
cleavagePositionNumberDF = pd.DataFrame(cleavagePositionNumber)
combinedWithPositionNumber = pd.DataFrame.join(combinedDF1, cleavagePositionNumberDF)

#Determine if protein is Lipoprotein
#Check if Predition has the lipoprotein in the results
isLipoproteinList = []

for eachOrganism in predictionSignalPeptide:
    findLipoprotein = eachOrganism.find('Lipoprotein')
    if findLipoprotein != -1:
        isLipoproteinList.append('Yes')
    else:
        isLipoproteinList.append('No')

#Create DateFrame on whether Organism is a lipoprotein
isLipoproteinResults = {'Is Lipoprotein': isLipoproteinList}
isLipoproteinResultsDF = pd.DataFrame(isLipoproteinResults)
combinedWithIsLipoprotein = pd.DataFrame.join(combinedWithPositionNumber, isLipoproteinResultsDF)

combineWithUniProtID = pd.DataFrame.join(combinedWithIsLipoprotein, organismUniProtIDDF)

#Keys for DataFrames
#Sort DataFrame by UniProtID in ascending order
signalP_Keys = combineWithUniProtID.keys()
keySignalP = signalP_Keys[-1] #Last element is UniProtID
combineWithUniProtID.sort_values(by=keySignalP, inplace=True)

#print(combineWithUniProtID[['Full_Sequence_Fasta']])
#print(combinedFastaDF3[['Full_Sequence_Fasta']])

#Merge Dataframes with matching UniProtID columns as Key
#combineWithUniProtID and combinedWithFastaDF2
finalDF = combinedFastaDF3.merge(combineWithUniProtID, left_on=keyFasta, right_on=keySignalP)

#Add Signal Peptide Sequence
finalDF_Keys = finalDF.keys()
cleavagePositionTemp = finalDF[['Position']]
fullSequenceTemp = finalDF[['Full_Sequence_Fasta']]

sequenceAndCleavage = pd.DataFrame.join(fullSequenceTemp, cleavagePositionTemp)

signalPeptideSequence = []
for eachOrganism in sequenceAndCleavage.itertuples():
    site = int(eachOrganism[-1])
    signalPeptideSequence.append(eachOrganism[1][0:site])

signalPeptideSeqDF = pd.DataFrame(signalPeptideSequence)
addedSignalPeptideDF = pd.DataFrame.join(finalDF, signalPeptideSeqDF)

#Export to Excel
addedSignalPeptideDF.to_excel(r'/Users/jeffreylai/Developer/Python/UniprotDownload/' + exportedExcelFileName)

#If message is seen in console, then the program completed
print("Process Completed ; SignalP JSON Results Parsed ; Check for Excel file: " + exportedExcelFileName)
