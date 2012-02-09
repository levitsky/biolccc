from pyteomics import biolccc

peptide = 'Ac-PEPTIDE-NH2'

averageMass = biolccc.calculateAverageMass(
    peptide, biolccc.rpAcnFaRod)
monoisotopicMass = biolccc.calculateMonoisotopicMass(
    peptide, biolccc.rpAcnFaRod)

print 'The average mass of', peptide, 'is', averageMass, 'Da'
print 'The monoisotopic mass of', peptide, 'is', monoisotopicMass, 'Da'
