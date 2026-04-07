from pyteomics import biolccc
peptide = 'PEPTIDE'

parsedSequence = biolccc.parseSequence(peptide, biolccc.rpAcnTfaChain)

for chemicalGroup in parsedSequence:
    print(chemicalGroup.name())
