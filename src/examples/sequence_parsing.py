from pyteomics import biolccc
peptide = 'PEPTIDE'

parsedSequence = biolccc.parseSequence(peptide)

for chemicalGroup in parsedSequence:
    print chemicalGroup.name()
