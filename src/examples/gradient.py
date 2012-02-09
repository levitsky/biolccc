from pyteomics import biolccc

peptide = 'Ac-PEPTIDE-NH2'
myChromoConditions = biolccc.ChromoConditions()
myGradient = biolccc.Gradient()
myGradient.addPoint(0.0, 5.0)
myGradient.addPoint(20.0, 5.0)
myGradient.addPoint(60.0, 45.0)
myGradient.addPoint(65.0, 100.0)
myGradient.addPoint(85.0, 100.0)
myChromoConditions.setGradient(myGradient)

RT = biolccc.calculateRT(peptide,
         biolccc.rpAcnFaRod,
         myChromoConditions)
print 'The retention time of', peptide, 'in the custom gradient is',RT
