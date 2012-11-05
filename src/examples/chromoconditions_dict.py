from pyteomics import biolccc

myChromoConditions = biolccc.ChromoConditions()
print myChromoConditions.keys()

myChromoConditions['columnLength'] = 100.0
myChromoConditions['columnDiameter'] = 0.1
myChromoConditions['columnPoreSize'] = 300.0
myChromoConditions['secondSolventConcentrationA'] = 5.0
myChromoConditions['secondSolventConcentrationB'] = 80.0
myChromoConditions['gradient'] = biolccc.Gradient(0.0, 90.0, 60.0)
myChromoConditions['flowRate'] = 0.0005

peptide = 'Ac-PEPTIDE-NH2'
RT = biolccc.calculateRT(peptide,
    biolccc.rpAcnFaRod,
    myChromoConditions)
print 'The retention time of', peptide, 'is', RT
