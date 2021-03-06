from pyteomics import biolccc 

for label, chemicalGroup in biolccc.rpAcnFaRod['chemicalGroups'].items():
    print 'Name', chemicalGroup['name']
    print 'Label', chemicalGroup['label']
    print 'Bind energy', chemicalGroup['bindEnergy']
    print 'Average mass', chemicalGroup['averageMass']
    print 'Monoisotopic mass', chemicalGroup['monoisotopicMass']
    print ''

print 'More simple syntax:'
print biolccc.rpAcnFaRod
