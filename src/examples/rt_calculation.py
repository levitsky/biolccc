from pyteomics import biolccc

peptide = 'Ac-PEPTIDE-NH2'
RT = biolccc.calculateRT(peptide,
    biolccc.rpAcnFaRod,
    biolccc.standardChromoConditions)
print 'The retention time of', peptide, 'is', RT
