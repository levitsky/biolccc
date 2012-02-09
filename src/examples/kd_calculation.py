from pyteomics import biolccc

peptide = 'Ac-PEPTIDE-NH2'

kd = biolccc.calculateKd(
    peptide, # the peptide sequence 
    15.0,    # the concentration of the second solvent, %
    biolccc.rpAcnFaRod, # the chemical basis
    100.0)   # the size of the pores, angstroms

print 'The coefficient of distribution of', peptide, 'is', kd
