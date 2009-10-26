#include <iostream>

#include "peptidemethods.h"
#include "math.h"

#define MEMORY_ERROR   -2.0
#define PARSING_ERROR  -3.0
#define PORESIZE_ERROR -4.0
#define GRADIENT_ERROR -5.0

double PeptideMethods :: calculateRTBioLCCC(const std::string &sequence,
                                            const ChemicalBasis &chemBasis,
                                            const ChromoConditions &conditions
) {
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;    
    Terminus CTerminus;
    std::vector<double> peptideEnergyProfile;
    
   if ( parseSequence(sequence, 
                    chemBasis,
                    &parsedPeptideStructure,
                    &NTerminus,
                    &CTerminus,
                    &peptideEnergyProfile) )
    {
        return calculateRTBioLCCC (peptideEnergyProfile,
                                chemBasis,
                                conditions);
    }
    else {
        return PARSING_ERROR;
    }
}

double PeptideMethods :: calculateRTBioLCCC(const std::vector<double> &peptideEnergyProfile,
                                            const ChemicalBasis &chemBasis,
                                            const ChromoConditions &conditions
) {
        
    // Calculating column volumes 
    double volumeLiquidPhase = conditions.columnDiameter() * conditions.columnDiameter() * 3.1415 * conditions.columnLength() / 4.0 / 1000.0 * (conditions.columnPorosity()-conditions.columnVpToVtot());
    double volumePore = conditions.columnDiameter() * conditions.columnDiameter() * 3.1415 * conditions.columnLength() / 4.0 / 1000.0 * conditions.columnVpToVtot();    
    
    // Recalculating dV. By default dV is calculated as flow rate divided by 20.
    double dV;
    if ( conditions.dV()!= 0 ) {
        dV = conditions.dV();
    }
    else {
        dV = conditions.flowRate() / 20.0;
    }

    // Because of the features of the BioLCCC model, size of a pore should be more than 20A.
    if (conditions.columnPoreSize() <=15) {
        return PORESIZE_ERROR;
    }

    // A gradient should start at time 0.
    if ( conditions.beginGradient()->first != 0 ) {
        return GRADIENT_ERROR;
    }
    
    // A gradient should contain at least two points.
    if ( ( conditions.beginGradient() == conditions.endGradient() ) || ( conditions.beginGradient() == --conditions.endGradient() ) ) {
        return GRADIENT_ERROR;
    }
    
    // Converting an x-coordinate of gradient to the dimension of iterations.
    std::vector<std::pair<int, double> > convertedGradient;
    double previousTimePoint = conditions.beginGradient()->first;
    
    for (std::vector<std::pair<double, double> >::const_iterator currentGradientPoint = conditions.beginGradient();
        currentGradientPoint < conditions.endGradient();
        currentGradientPoint++
    ) {
        // Gradient should be ascending in a time.
        if ( (previousTimePoint != conditions.beginGradient()->first ) && (previousTimePoint >= currentGradientPoint->first ) ) {
            return GRADIENT_ERROR;
        }
        convertedGradient.push_back(std::pair<int,double>(currentGradientPoint -> first * conditions.flowRate() / dV,
                                                          currentGradientPoint -> second));
        
        previousTimePoint = currentGradientPoint->first;
    }
    
    double S = 0.0; // A part of column passed by molecules. When it exceeds 1.0, molecule elute from a column.
    int j = 0;      // A current iteration number. 
    double secondSolventConcentration = 0.0;


    std::vector<std::pair<int, double> >::const_iterator currentGradientPoint  = convertedGradient.begin();
    std::vector<std::pair<int, double> >::const_iterator previousGradientPoint = convertedGradient.begin();
    while ( S < 1.0 ) {
        j++;
        if ( j > currentGradientPoint->first ) {
            // If j exceeds the last point of a gradient, the value of a second solvent concentration is calculated by a prolongation of a last gradient section.
            if ( currentGradientPoint != --convertedGradient.end() ) {
                previousGradientPoint = currentGradientPoint;
                ++currentGradientPoint;
                
                // We could calculate isocratic part of gradient more efficiently due to constant Kd.
                if (currentGradientPoint->second == previousGradientPoint -> second) {
                    // One case is that a peptide elutes during this section or the section is the last.
                    if ( ( ( 1 - S ) / dV * calculateKdBioLCCC (peptideEnergyProfile, chemBasis, conditions.secondSolvent(), secondSolventConcentration, conditions.columnPoreSize(), conditions.calibrationParameter(),conditions.temperature()) * volumePore < ( currentGradientPoint->first - j + 1 ) ) ||
                    (currentGradientPoint == --convertedGradient.end() ) ) {
                         j += ceil( ( 1 - S ) / dV * calculateKdBioLCCC (peptideEnergyProfile, chemBasis, conditions.secondSolvent(), secondSolventConcentration, conditions.columnPoreSize(), conditions.calibrationParameter(),conditions.temperature())  * volumePore ) - 1;
                         break;
                    }
                    // Another case is that this section is not enough for a peptide to elute.
                    else {
                        S += dV / calculateKdBioLCCC (peptideEnergyProfile, chemBasis, conditions.secondSolvent(), secondSolventConcentration, conditions.columnPoreSize(), conditions.calibrationParameter(),conditions.temperature()) / volumePore * ( currentGradientPoint->first - j + 1 );
                        j = currentGradientPoint->first;
                    }
                }
            }
        }
        secondSolventConcentration = currentGradientPoint->second - (currentGradientPoint->second - previousGradientPoint->second) / (currentGradientPoint->first - previousGradientPoint->first ) * ( currentGradientPoint->first - j );
        //std::cout << j << " " << secondSolventConcentration << "<br>";
        S += dV / calculateKdBioLCCC (peptideEnergyProfile, chemBasis, conditions.secondSolvent(), secondSolventConcentration, conditions.columnPoreSize(), conditions.calibrationParameter(),conditions.temperature()) / volumePore;
        //std::cout << secondSolventConcentration << " " << 1 / (1 + calculateKdBioLCCC (peptideEnergyProfile, chemBasis, conditions.secondSolvent(), secondSolventConcentration, conditions.columnPoreSize(), conditions.calibrationParameter(),conditions.temperature()) ) << "\n" << "<br>";
    }

    double RT = 0;
    
    
    
    RT += j * dV / conditions.flowRate() + conditions.delayTime() + volumeLiquidPhase / conditions.flowRate();
    //RT = j;
    return RT;
}

double PeptideMethods :: calculateKdBioLCCC (const std::string &sequence,
                                            const ChemicalBasis &chemBasis,
                                            const std::string &secondSolvent,
                                            const double secondSolventConcentration,
                                            const double columnPoreSize,
                                            const double calibrationParameter,
                                            const double temperature
        
) {
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;    
    Terminus CTerminus;
    std::vector<double> peptideEnergyProfile;
    
   if ( parseSequence(sequence, 
                    chemBasis,
                    &parsedPeptideStructure,
                    &NTerminus,
                    &CTerminus,
                    &peptideEnergyProfile) )
    {
        return calculateKdBioLCCC (peptideEnergyProfile,
                                chemBasis,
                                secondSolvent,
                                secondSolventConcentration,
                                columnPoreSize,
                                calibrationParameter,
                                temperature);
    }
    else {
        return PARSING_ERROR;
    }
}                                

double PeptideMethods :: calculateKdBioLCCC (const std::vector<double> &peptideEnergyProfile,
                                            const ChemicalBasis &chemBasis,
                                            const std::string secondSolvent,
                                            const double secondSolventConcentration,
                                            const double columnPoreSize,
                                            const double calibrationParameter,
                                            const double temperature
) {
    // At first, we need to convert energy profile to a profile of distribution probabilities. Probability = exp (E_effective),
    // where E_effective = E_of_residue - Eab,
    // and Eab is an energy of binding for a binary solvent and Eab = Ea + ln ( 1 + Nb + Nb*exp (Ea - Eb) )
    // also corrections of energies due to temperature (energies in exponents are scaled to the new temperature) 
    // and column aging (calibration parameter) are introduced.
    
    // Due to preliminary scaling the binding energy of water equals zero.
    double Q = exp ( ( 0 + chemBasis.secondSolventBindEnergy() ) * 293.0 / temperature );
    double Nb = 0;
    double Eab = 0;
   
    // A Boltzmann factor is an exponent of an energy of interaction between an amino acid residue and a solid phase divided by a temperature * Boltzmann's constant. An energy unit is a Boltzmann's constant * 293.0 Kelvins.
    // This probability is used later in transition matrix.
    std::vector<double> boltzmannFactorProfile;
    for ( std::vector<double>::const_iterator residueEnergy = peptideEnergyProfile.begin();
        residueEnergy != peptideEnergyProfile.end();
        residueEnergy ++
    ) {
        //std::cout << *residueEnergy << " ";
        Nb = secondSolventConcentration * 1.91 / ( secondSolventConcentration * 1.91 + (100.0 - secondSolventConcentration) * 5.56 );
        Eab = 0.0 + 1.0 / calibrationParameter * log( 1.0 - Nb + Nb * Q );
        boltzmannFactorProfile.push_back( exp ( calibrationParameter * ( *residueEnergy - Eab ) * 293.0 / temperature ) );
    }
    
    // The density vector correspond to a probability of n-th residue to be in a certain layer between pore walls.
    // The transition matrix used to calculate a density vector of n-th particle from a density vector of (n-1)-th particle.
    // The density buffer vector is used during matrix calculations.
    double *density;
    double *transitionMatrix;
    double *densityBuffer;
    
    // PoreSteps is a number of nodes in a lattice between two walls. Because of the features of a following calculation it should be more than 2.
    const int poreSteps = (int) columnPoreSize / 10.0 + 1;
    if (poreSteps <=2) {
        return PORESIZE_ERROR;
    }
    
    // Memory managment.
    try {
        density = new double[poreSteps];
        densityBuffer = new double[poreSteps];
        transitionMatrix = new double[poreSteps*poreSteps];
    }
    catch (...) {
        return MEMORY_ERROR;
    }

    // Now we need to construct a density vector for the last aminoacid residue.
    // A density is distributed uniformely over all layers of the lattice, except for the wall layers. There a density is multiplied by Boltzmann factor due to interaction of a residue with a solid phase.
    density[0] = boltzmannFactorProfile[ boltzmannFactorProfile.size() - 1 ];
    for (int i = 1; i < poreSteps - 1 ; i++) {
        density[i] = 1;
    }
    density[ poreSteps - 1 ] = boltzmannFactorProfile[ boltzmannFactorProfile.size() - 1 ];

    // Than we construct a basis for the transition matrix. The basis is a diagonal matrix with 4.0/6.0 on a main diagonal and 1.0/6.0 on side diagonals.
    
    // Filling the first row.
    transitionMatrix[0] = 4.0/6.0;
    transitionMatrix[1] = 1.0/6.0;
    for (int i = 2; i < poreSteps; i++) {
        transitionMatrix[i] = 0.0;
    }

    // Filling from 2nd to (n-1)th rows.
    for (int i = 1; i < poreSteps - 1; i++) {
        for (int j = 0; j < poreSteps; j++) {
            switch ( j - i + 1 ) {
                case 0: {
                    transitionMatrix[j + poreSteps * i] = 1.0/6.0;
                    break;
                }
                case 1: {
                    transitionMatrix[j + poreSteps * i] = 4.0/6.0;
                    break;
                }
                case 2: {
                    transitionMatrix[j + poreSteps * i] = 1.0/6.0;
                    break;
                }
                default:
                    transitionMatrix[j + poreSteps * i] = 0.0;
            }
        }
    }
    
    // Filling the n-th row.
    for (int i = poreSteps * (poreSteps - 1); i < ( poreSteps * poreSteps - 2); i++) {
        transitionMatrix[i] = 0.0;
    }
    transitionMatrix[poreSteps * poreSteps - 2] = 1.0/6.0;
    transitionMatrix[poreSteps * poreSteps - 1] = 4.0/6.0;
    
    // On the each step we calculate a density vector for the n-th amino acid residue by a multiplication of a transition matrix and a density vector for the (n-1)th residue.
    for ( std::vector<double>::const_reverse_iterator residueBoltzmannFactor = ++boltzmannFactorProfile.rbegin();
        residueBoltzmannFactor != boltzmannFactorProfile.rend();
        residueBoltzmannFactor++
    ) {
        // Elements of the first and the last rows of a transition matrix are modified by Boltzmann factor.
        transitionMatrix[0] = 4.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[1] = 1.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[poreSteps*poreSteps - 1] = 4.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[poreSteps*poreSteps - 2] = 1.0 / 6.0 * (*residueBoltzmannFactor);

        // Zeroing the calculation buffer.
        for (int i = 0; i < poreSteps; i++) {
            densityBuffer[i] = 0.0;
        }

        // Multiplying a transition matrix by a density vector. Result is stored in a buffer vector.
        for (int i = 0; i < poreSteps; i++) {
            for (int j = 0; j < poreSteps; j++) {
                densityBuffer[i] = densityBuffer[i] + density[j] * transitionMatrix[ j + i * poreSteps];
            }
        } 

        // Transferring results in a density vector.
        for (int i = 0; i < poreSteps; i++)  {
            density[i] = densityBuffer[i];
        }
    }

    // Finally, Kd is calculated as a sum of elements of a density vector. It's normalized on a size of a lattice.
    double Kd=0;
    for (int k=0;k<poreSteps;k++) {
        Kd += density[k];
    }
    Kd = Kd / poreSteps;

    // Cleaning memory.
    try {
        delete density;
        delete densityBuffer;
        delete transitionMatrix;
    }
    catch (...) {
        return MEMORY_ERROR;
    }

    return Kd;
}                                            

double PeptideMethods :: calculateAverageMass(const std::string &sequence,
                                            const ChemicalBasis &chemBasis
){
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;
    Terminus CTerminus;
    double peptideAverageMass = 0;
    
    if ( parseSequence (sequence, chemBasis, &parsedPeptideStructure, &NTerminus, &CTerminus) ) {
        for ( std::vector<Aminoacid>::const_iterator i = parsedPeptideStructure.begin(); i < parsedPeptideStructure.end(); i++) {
            peptideAverageMass += i -> averageMass();
        }
        peptideAverageMass += NTerminus.averageMass();
        peptideAverageMass += CTerminus.averageMass();
    }
    
    return peptideAverageMass;
}
                                        
double PeptideMethods :: calculateMonoisotopicMass(const std::string &sequence,
                                                const ChemicalBasis &chemBasis
) {
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;
    Terminus CTerminus;
    double monoisotopicMass = 0;
    
    if ( parseSequence (sequence, chemBasis, &parsedPeptideStructure, &NTerminus, &CTerminus) ) {
        for ( std::vector<Aminoacid>::const_iterator i = parsedPeptideStructure.begin(); i < parsedPeptideStructure.end(); i++) {
            monoisotopicMass += i -> monoisotopicMass();
        }
        monoisotopicMass += NTerminus.monoisotopicMass();
        monoisotopicMass += CTerminus.monoisotopicMass();
    }
    
    return monoisotopicMass;
}

bool PeptideMethods :: parseSequence(const std::string &source, 
                                    const ChemicalBasis &chemBasis,
                                    std::vector<Aminoacid> *parsedPeptideStructure,
                                    Terminus *NTerminus,
                                    Terminus *CTerminus,
                                    std::vector<double> *peptideEnergyProfile
) {
    parsedPeptideStructure -> clear();
    
    // At first we need to strip the sequence from adjacent aminoacids.
    // If a source sequence contains them, it should contain two dots, so we need the part of sequence between them.
    std::size_t firstDotPosition = 0;
    std::size_t secondDotPosition = 0;
    
    // We'll use this variable to contain peptide sequence without adjacent aminoacids. 
    std::string strippedSource = source;
    
    firstDotPosition = source.find(".");
    
    if ( firstDotPosition != std::string::npos ) {
        secondDotPosition = source.find(".", firstDotPosition+1);
        if ( secondDotPosition != std::string::npos ) {  
      
            // If a source sequence contains more that two dots, it's broken.
            if ( source.find(".", secondDotPosition+1) != std::string::npos ) {
                return false;
            }
            else {
                strippedSource = source.substr(firstDotPosition+1, secondDotPosition - firstDotPosition - 1);
            }
        }
        // If a source sequence contains only one dot, it's broken.
        else {
            return false;
        }
    }
    
    
    // Than goes parsing.
    std::size_t NTerminusPosition = 0;
    
    // First we need to check the stripped source sequence for a N-Terminal group.
    *NTerminus = chemBasis.defaultNTerminus()->second ;
    for (std::map<std::string,Terminus>::const_iterator NTermIterator = chemBasis.beginNTerminus(); 
        NTermIterator != chemBasis.endNTerminus(); 
        NTermIterator++
    ) {
        if ( strippedSource.find( NTermIterator->second.label() ) != std::string::npos ) {
            *NTerminus = NTermIterator->second;
            NTerminusPosition = NTermIterator->second.label().size();
        }
    }
    
    // Then we need to found the location of the C-Terminus.
    *CTerminus = chemBasis.defaultCTerminus()->second;
    std::size_t CTerminusPosition;
    CTerminusPosition = strippedSource.find("-", NTerminusPosition);
    if ( CTerminusPosition != std::string::npos ){
        
        // The sequence should not contain hyphens after C-terminal group.
        if (strippedSource.find("-", CTerminusPosition+1) != std::string::npos) {
            return false;
        }
        
        // Searching for known C-terminal groups.
        for (std::map<std::string,Terminus>::const_iterator CTermIterator = chemBasis.beginCTerminus(); 
            CTermIterator != chemBasis.endCTerminus(); 
            CTermIterator++
        ) {
            if ( strippedSource.find( CTermIterator->second.label(), CTerminusPosition ) != std::string::npos ) {
                *CTerminus = CTermIterator->second;
            }
        }
    }
    else {
        CTerminusPosition = strippedSource.size();
    }
        
    // At this step we obtain the sequence of a peptide without adjacent aminoacids and terminal groups.
    strippedSource = strippedSource.substr(NTerminusPosition, CTerminusPosition-NTerminusPosition);


    // We need to check whether it contains any non-letter characters.
    for (std::size_t i=0; i<strippedSource.size(); i++) {
        if ( !( ( ( int(strippedSource[i]) >= int('a') ) && ( int(strippedSource[i]) <= int('z') ) ) ||
              ( ( int(strippedSource[i]) >= int('A') ) && ( int(strippedSource[i]) <= int('Z') ) ) ) ) {
            return false;
        }
    }

    // Then we divide the whole sequence on aminoacids. 
    bool aminoacidFound;
    size_t curPos = 0;
    while (curPos < strippedSource.size()) {
        aminoacidFound = false;
        for (std::map<std::string,Aminoacid>::const_iterator currentAminoacid = chemBasis.beginAminoacid(); 
            currentAminoacid != chemBasis.endAminoacid(); 
            currentAminoacid++
        ) {
            if ( strippedSource.compare(curPos, currentAminoacid->second.label().size(), currentAminoacid->second.label() ) == 0 ) {
                curPos += currentAminoacid->second.label().size();
                //std::cout << currentAminoacid->second.name() << "\n";
                parsedPeptideStructure -> push_back(currentAminoacid->second);
                aminoacidFound = true;
                break;
            }
        }
  
        if (!aminoacidFound) {
            return false;
        }
    }
    
    // Finally, we build an energy profile if it was defined.
    if (peptideEnergyProfile != NULL)
    {
        peptideEnergyProfile -> clear();
        for (std::vector<Aminoacid>::const_iterator currentAminoacid = parsedPeptideStructure -> begin();
            currentAminoacid != parsedPeptideStructure -> end();
            currentAminoacid++
        ) {
            peptideEnergyProfile -> push_back(currentAminoacid->bindEnergy());
        }
        
        // Modifing energies of terminal aminoacid residues.
        *(peptideEnergyProfile->begin()) = *(peptideEnergyProfile->begin()) + NTerminus->bindEnergy();
        *(--peptideEnergyProfile->end()) = *(--peptideEnergyProfile->end()) + CTerminus->bindEnergy();
    }

    return true;
}

bool PeptideMethods :: calculatePeptideProperties(const std::string &sequence,
                                                        const ChemicalBasis &chemBasis,
                                                        const ChromoConditions &conditions,
                                                        double *RTBioLCCC,
                                                        double *averageMass,
                                                        double *monoisotopicMass
) {
   std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;    
    Terminus CTerminus;
    std::vector<double> peptideEnergyProfile;
    
   if ( parseSequence(sequence, 
                    chemBasis,
                    &parsedPeptideStructure,
                    &NTerminus,
                    &CTerminus,
                    &peptideEnergyProfile) )
    {
        *RTBioLCCC = calculateRTBioLCCC (peptideEnergyProfile, chemBasis, conditions);
        
        *averageMass = 0;
        for ( std::vector<Aminoacid>::const_iterator i = parsedPeptideStructure.begin(); i < parsedPeptideStructure.end(); i++) {
            *averageMass += i -> averageMass();
        }
        *averageMass += NTerminus.averageMass();
        *averageMass += CTerminus.averageMass();
        
        *monoisotopicMass = 0;
        for ( std::vector<Aminoacid>::const_iterator i = parsedPeptideStructure.begin(); i < parsedPeptideStructure.end(); i++) {
            *monoisotopicMass += i -> monoisotopicMass();
        }
        *monoisotopicMass += NTerminus.monoisotopicMass();
        *monoisotopicMass += CTerminus.monoisotopicMass();
        
        return true;
    }
    else {
        return false;
    }
}

