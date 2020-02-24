class UniProt:
    urlBase = "https://www.uniprot.org/uniprot/"
    urlFileType = ".fasta"

    def __init__(self, organism, locusTag, uniprotID, isLipoprotein):
        self.organism = organism
        self.locusTag = locusTag
        self.uniprotID = uniprotID
        self.isLipoprotein = isLipoprotein
