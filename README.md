# UniProtPython
#**Python Program**
##Python Program to extract data from UniProt FASTA File and SignalP Server

SignalP Server:
http://www.cbs.dtu.dk/services/SignalP/

UniProt:
https://www.uniprot.org/uniprot/

Files:

Version 1: (3/6/20)

 - UniProtObject.py :
  Data Model
 - UniProtSearch.py :
  Obtain results from UniProt Protein FASTA File
 - SignalPjsonDecoder.py :
  Analyze results from SignalP 5.0 Server ; Input is JSON File

Version 2: (5/1/20)

 - UniprotObjectFasta.py :
  Data Model
 - UniProtSearchVer2.py :
  Example to parse UniProt Protein FASTA File , includes HTTP request with Uniprot (Optional Component)
 - SignalPJsonDecoder2.py :
  Extract data from UniProt Protein FASTA file and SignalP (Version 5) JSON Results and exports results into Excel Spreadsheet


Version 3: (6/8/20)

 - UniRef90Object.py:
  Data Model for New UniRef90 fasta format  
 - SignalPJsonDecoder3.py :
  Extract data from UniRef90 Protein FASTA file and SignalP (Version 5) JSON results and expors results into Excel Spreadsheet
