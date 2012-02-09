from pyteomics import biolccc

print biolccc.calculateRT('QWERTYIPASDFGHKLCVNM', biolccc.rpAcnTfaChain,
    biolccc.standardChromoConditions)

# Using 21 interpolating points.
print biolccc.calculateRT('QWERTYIPASDFGHKLCVNM', biolccc.rpAcnTfaChain,
    biolccc.standardChromoConditions, 21)
