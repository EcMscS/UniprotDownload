#Modules needed for Python program to work
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
#-------------------------------------------------------------

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

for each in uniprot_objects:
    info = each.organism + " " + each.uniprotID + " " + each.description + " " + each.fullSequence
#---------------------------------------

with open(path + fileName, "r") as read_file:
    data = json.load(read_file)

#Prints all sequences from JSON from top level
print(data['SEQUENCES'])

#Parse JSON Function
organism = extract_values(data, 'Name')
cleavageSite = extract_values(data, 'CS_pos')
predictionSignalPeptide = extract_values(data, 'Prediction')

#Create 2 DataFrames
organismList = {'organism': organism}
cleavageSiteList = {'cleavage_site': cleavageSite}

organismDF = pd.DataFrame(organismList)
cleavageSiteDF = pd.DataFrame(cleavageSiteList)
predictionListDF = pd.DataFrame(predictionSignalPeptide)

combinedDF1 = pd.DataFrame.join(organismDF, cleavageSiteDF)
combinedFinalDF = pd.DataFrame.join(combinedDF1, predictionListDF)

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

#Add fullSequence to dataFrame
fullSequenceOfOrganism = {'Full_Sequence':extractedProteinSequence}
fullSequenceOfOrganismDF = pd.DataFrame(fullSequenceOfOrganism)
addedFullSequence = pd.DataFrame.join(combinedWithPositionNumber, fullSequenceOfOrganismDF)

#Cut sequences according to cleavage site
#Use extractedLine2FromFasta to create new signal peptide sequence
signalPeptideSequence = []

#PROBLEM HERE
# 3 is position and 4 is Full Sequence
for eachOrganism in addedFullSequence.itertuples():
    site = int(eachOrganism[3])
    signalPeptideSequence.append(eachOrganism[4][0:site])

#Add signal peptide sequence to dataFrame
signalPeptideSequenceOfOrganism = {'Signal_Sequence':signalPeptideSequence}
signalPeptideSequenceDF = pd.DataFrame(signalPeptideSequenceOfOrganism)
addedSignalPeptideSequence = pd.DataFrame.join(addedFullSequence, signalPeptideSequenceDF)

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
combinedWithIsLipoprotein = pd.DataFrame.join(addedSignalPeptideSequence, isLipoproteinResultsDF)

#Export to Excel
combinedWithIsLipoprotein.to_excel(r'/Users/jeffreylai/Developer/Python/UniprotDownload/' + exportedExcelFileName)


#If message is seen in console, then the program completed
print("Process Completed ; SignalP JSON Results Parsed ; Check for Excel file: " + exportedExcelFileName)
