class UniProt:
    urlBase = "https://www.uniprot.org/uniprot/"
    urlFileType = ".fasta"

    def __init__(self, organism, uniprotID, description, fullSequence):
        self.organism = organism
        self.uniprotID = uniprotID
        self.description = description
        self.fullSequence = fullSequence
