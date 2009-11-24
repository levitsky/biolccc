#include <iostream>

#include "boost/foreach.hpp"

#include "BioLCCC.h"
#include "math.h"

#define MEMORY_ERROR   -2.0
#define PARSING_ERROR  -3.0
#define PORESIZE_ERROR -4.0
#define GRADIENT_ERROR -5.0

namespace BioLCCC {

namespace {

// Auxiliary functions that shouldn't be exposed to user at this
// point.
bool parseSequence(
    const std::string &source, 
    const ChemicalBasis &chemBasis,
    std::vector<Aminoacid> *parsedPeptideStructure,
    Terminus *NTerminus,
    Terminus *CTerminus,
    std::vector<double> *peptideEnergyProfile
) {
    parsedPeptideStructure -> clear();
    
    // At first we need to strip the sequence from adjacent aminoacids.
    // If a source sequence contains them, it should contain two dots, so we 
    // need the part of sequence between them.
    std::size_t firstDotPosition = 0;
    std::size_t secondDotPosition = 0;
    
    // We'll use this variable to contain peptide sequence without adjacent 
    // aminoacids. 
    std::string strippedSource = source;
    
    firstDotPosition = source.find(".");
    
    if (firstDotPosition != std::string::npos) {
        secondDotPosition = source.find(".", firstDotPosition+1);
        if (secondDotPosition != std::string::npos) {  
      
            // If a source sequence contains more that two dots, it's broken.
            if (source.find(".", secondDotPosition+1) != std::string::npos) {
                return false;
            }
            else {
                strippedSource = source.substr(firstDotPosition+1, 
                    secondDotPosition - firstDotPosition - 1);
            }
        }
        // If a source sequence contains only one dot, it's broken.
        else {
            return false;
        }
    }
    
    
    // Than goes parsing.
    std::size_t NTerminusPosition = 0;
    
    // First we need to check the stripped source sequence for 
    // the N-Terminal group.
    *NTerminus = chemBasis.defaultNTerminus();
    std::pair<std::string,Terminus> NTerminusIterator;
    BOOST_FOREACH(NTerminusIterator, chemBasis.NTermini()) {
        if (strippedSource.find(NTerminusIterator.second.label()) != 
                std::string::npos) {
            *NTerminus = NTerminusIterator.second;
            NTerminusPosition = NTerminusIterator.second.label().size();
        }
    }
    
    // Then we need to found the location of the C-Terminus.
    *CTerminus = chemBasis.defaultCTerminus();
    std::size_t CTerminusPosition;
    CTerminusPosition = strippedSource.find("-", NTerminusPosition);
    if (CTerminusPosition != std::string::npos){
        
        // The sequence should not contain hyphens after C-terminal group.
        if (strippedSource.find("-", CTerminusPosition+1) != std::string::npos){
            return false;
        }
        
        // Searching for known C-terminal groups.
        std::pair<std::string,Terminus> CTerminusIterator;
        BOOST_FOREACH(CTerminusIterator, chemBasis.CTermini()) {
            if (strippedSource.find(CTerminusIterator.second.label(), 
                CTerminusPosition ) != std::string::npos) {
                *CTerminus = CTerminusIterator.second;
            }
        }
    }
    else {
        CTerminusPosition = strippedSource.size();
    }
        
    // At this step we obtain the sequence of a peptide without adjacent 
    // aminoacids and terminal groups.
    strippedSource = strippedSource.substr(NTerminusPosition, 
                     CTerminusPosition-NTerminusPosition);

    // We need to check whether it contains any non-letter characters.
    for (std::size_t i=0; i<strippedSource.size(); i++) {
        if (!(((int(strippedSource[i]) >= int('a')) && 
               (int(strippedSource[i]) <= int('z'))) ||
              ((int(strippedSource[i]) >= int('A')) && 
               (int(strippedSource[i]) <= int('Z'))))) {
            return false;
        }
    }

    // Then we divide the whole sequence on aminoacids. 
    bool aminoacidFound;
    size_t curPos = 0;
    while (curPos < strippedSource.size()) {
        aminoacidFound = false;
        std::pair<std::string,Aminoacid> aminoacidIterator;
        BOOST_FOREACH(aminoacidIterator, chemBasis.aminoacids()) {
            if (strippedSource.compare(curPos, 
                    aminoacidIterator.second.label().size(), 
                    aminoacidIterator.second.label()) == 0) {
                curPos += aminoacidIterator.second.label().size();
                parsedPeptideStructure -> push_back(aminoacidIterator.second);
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
        peptideEnergyProfile->clear();
        BOOST_FOREACH(Aminoacid currentAminoacid, *parsedPeptideStructure){
            peptideEnergyProfile->push_back(currentAminoacid.bindEnergy());
        }
        
        // Modifing energies of terminal aminoacid residues.
        *(peptideEnergyProfile->begin()) = *(peptideEnergyProfile->begin()) +
            NTerminus->bindEnergy();
        *(--peptideEnergyProfile->end()) = *(--peptideEnergyProfile->end()) + 
            CTerminus->bindEnergy();
    }

    return true;
}

double calculateKd(
    const std::vector<double> &peptideEnergyProfile,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double calibrationParameter,
    const double temperature
) {
    // At first, we need to convert energy profile to a profile of distribution
    // probabilities. Probability = exp(E_effective),
    // where E_effective = E_of_residue - Eab,
    // and Eab is an energy of binding for a binary solvent 
    // and Eab = Ea + ln ( 1 + Nb + Nb*exp (Ea - Eb) )
    // also corrections of energies due to temperature (energies in exponents 
    // are scaled to the new temperature) and column aging (calibration
    // parameter) are introduced.
    
    // Due to the preliminary scaling the binding energy of water equals zero.
    double Q = exp((0 + chemBasis.secondSolventBindEnergy()) * 
               293.0 / temperature);
    double Nb = 0;
    double Eab = 0;
   
    // A Boltzmann factor is an exponent of an energy of interaction between
    // an amino acid residue and a solid phase divided by a temperature * 
    // Boltzmann's constant. An energy unit is a Boltzmann's constant * 
    // 293.0 Kelvins. This probability is used later in the transition matrix.
    std::vector<double> boltzmannFactorProfile;
    for (std::vector<double>::const_iterator residueEnergy = 
             peptideEnergyProfile.begin();
         residueEnergy != peptideEnergyProfile.end();
         residueEnergy ++
    ) {
        //std::cout << *residueEnergy << " ";
        Nb = secondSolventConcentration * 1.91 / 
             (secondSolventConcentration * 1.91 + 
              (100.0 - secondSolventConcentration) * 5.56);
        Eab = 0.0 + 1.0 / calibrationParameter * log( 1.0 - Nb + Nb * Q );
        boltzmannFactorProfile.push_back(exp(calibrationParameter * 
            (*residueEnergy - Eab) * 293.0 / temperature));
    }
    
    // The density vector correspond to a probability of n-th residue to be in 
    // a certain layer between pore walls.
    // The transition matrix used to calculate a density vector of n-th
    // particle from a density vector of (n-1)-th particle.
    // The density buffer vector is used during matrix calculations.
    double *density;
    double *transitionMatrix;
    double *densityBuffer;
    
    // PoreSteps is a number of nodes in a lattice between two walls. Because of
    // the features of a following calculation it should be more than 2.
    const int poreSteps = (int) columnPoreSize / 10.0;
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
    // A density is distributed uniformely over all layers of the lattice, 
    // except for the wall layers. There a density is multiplied by Boltzmann 
    // factor due to interaction of a residue with a solid phase.
    density[0] = boltzmannFactorProfile[0];
    for (int i = 1; i < poreSteps - 1 ; i++) {
        density[i] = 1;
    }
    density[poreSteps - 1] = 
        boltzmannFactorProfile[0];

    // Than we construct a basis for the transition matrix. The basis is 
    // a diagonal matrix with 4.0/6.0 on the main diagonal and 1.0/6.0 on 
    // the side diagonals.
    
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
    for (int i = poreSteps * (poreSteps - 1); 
         i < (poreSteps * poreSteps - 2); 
         i++) {
        transitionMatrix[i] = 0.0;
    }
    transitionMatrix[poreSteps * poreSteps - 2] = 1.0/6.0;
    transitionMatrix[poreSteps * poreSteps - 1] = 4.0/6.0;
    
    // On the each step we calculate a density vector for the n-th aminoacid 
    // residue by multiplication of the transition matrix and the density vector
    // for the (n-1)th residue.
    for (std::vector<double>::const_iterator residueBoltzmannFactor =
             ++boltzmannFactorProfile.begin();
         residueBoltzmannFactor != boltzmannFactorProfile.end();
         residueBoltzmannFactor++
    ) {
        // Elements of the first and the last rows of the transition matrix are
        // modified by Boltzmann factor.
        transitionMatrix[0] = 4.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[1] = 1.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[poreSteps*poreSteps - 1] = 4.0 / 6.0 * 
            (*residueBoltzmannFactor);
        transitionMatrix[poreSteps*poreSteps - 2] = 1.0 / 6.0 * 
            (*residueBoltzmannFactor);

        // Zeroing the calculation buffer.
        for (int i = 0; i < poreSteps; i++) {
            densityBuffer[i] = 0.0;
        }

        // Multiplying the transition matrix by the density vector. The result 
        // is stored in the buffer vector.
        for (int i = 0; i < poreSteps; i++) {
            for (int j = 0; j < poreSteps; j++) {
                densityBuffer[i] = densityBuffer[i] + density[j] * 
                    transitionMatrix[j + i * poreSteps];
            }
        } 

        // Transferring results from the density vector.
        for (int i = 0; i < poreSteps; i++)  {
            density[i] = densityBuffer[i];
        }
    }

    // Finally, Kd is calculated as a sum of elements of the density vector. 
    // It is normalized to the size of the lattice.
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

double calculateRT(const std::vector<double> &peptideEnergyProfile,
    const ChromoConditions &conditions,
    const ChemicalBasis &chemBasis
) {
        
    // Calculating column volumes 
    double volumeLiquidPhase = conditions.columnDiameter() * 
        conditions.columnDiameter() * 3.1415 * conditions.columnLength() / 4.0 /
        1000.0 * (conditions.columnPorosity()-conditions.columnVpToVtot());
    double volumePore = conditions.columnDiameter() * 
        conditions.columnDiameter() * 3.1415 * conditions.columnLength() / 4.0 /
        1000.0 * conditions.columnVpToVtot();    
    
    // Recalculating dV. By default dV is calculated as flow rate divided by 20.
    double dV;
    if ( conditions.dV()!= 0.0 ) {
        dV = conditions.dV();
    }
    else {
        dV = conditions.flowRate() / 20.0;
    }

    // Because of the features of the BioLCCC model, size of the pore should be 
    // more than 20A.
    if (conditions.columnPoreSize() <= 15) {
        return PORESIZE_ERROR;
    }

    // A gradient should start at time 0.
    if (conditions.gradient().begin()->time() != 0.0) {
        return GRADIENT_ERROR;
    }
    
    // A gradient should contain at least two points.
    if (conditions.gradient().size() < 2) {
        return GRADIENT_ERROR;
    }
    
    // Converting the x-coordinate of the gradient to the scale of iterations
    // and the y-coordinate to the scale of second solvent concentration.
    std::vector<std::pair<int, double> > convertedGradient;
    
    BOOST_FOREACH(GradientPoint currentGradientPoint, conditions.gradient()) {
        convertedGradient.push_back(
                std::pair<int,double>(
                    floor(currentGradientPoint.time() * conditions.flowRate()
                        / dV),
                    (100.0 - currentGradientPoint.concentrationB()) / 100.0 *
                        conditions.secondSolventConcentrationA() + 
                        currentGradientPoint.concentrationB() / 100.0 *
                        conditions.secondSolventConcentrationB()));
    }
    
     // The part of a column passed by molecules. When it exceeds 1.0,
     // molecule elute from the column.
    double S = 0.0;
    // The current iteration number. 
    int j = 0;      
    double secondSolventConcentration = 0.0;


    std::vector<std::pair<int, double> >::const_iterator currentGradientPoint =
        convertedGradient.begin();
    std::vector<std::pair<int, double> >::const_iterator previousGradientPoint =
        convertedGradient.begin();
    while (S < 1.0) {
        j++;
        if (j > currentGradientPoint->first) {
            // If j exceeds the last point of a gradient, the value of the 
            // second solvent concentration is calculated by a prolongation of
            // the last gradient section.
            if (currentGradientPoint != --convertedGradient.end()) {
                previousGradientPoint = currentGradientPoint;
                ++currentGradientPoint;
                
                // We could calculate the isocratic part of gradient more 
                // efficiently due to the constancy of Kd.
                if (currentGradientPoint->second == 
                    previousGradientPoint->second) {
                    // One case is that a peptide elutes during this section or
                    // the section is the last.
                    if (((1.0 - S) / dV * 
                        calculateKd(peptideEnergyProfile, 
                            currentGradientPoint->second, 
                            chemBasis, 
                            conditions.columnPoreSize(), 
                            conditions.calibrationParameter(),
                            conditions.temperature()) *
                        volumePore < (currentGradientPoint->first - j + 1)) ||
                        
                        (currentGradientPoint == --convertedGradient.end()))
                    {
                        j += ceil((1.0 - S) / dV * 
                            calculateKd(peptideEnergyProfile, 
                                currentGradientPoint->second,
                                chemBasis,
                                conditions.columnPoreSize(),
                                conditions.calibrationParameter(),
                                conditions.temperature()) * volumePore);
                         break;
                    }
                    // Another case is that this section is not long enough for 
                    // a peptide to elute.
                    else {
                        S += dV / calculateKd(peptideEnergyProfile,
                                      currentGradientPoint->second, 
                                      chemBasis,
                                      conditions.columnPoreSize(), 
                                      conditions.calibrationParameter(),
                                      conditions.temperature()) /
                             volumePore * (currentGradientPoint->first -
                                           previousGradientPoint->first);
                        j = currentGradientPoint->first;
                    }
                }
            }
        }
        secondSolventConcentration = currentGradientPoint->second - 
            (currentGradientPoint->second - previousGradientPoint->second) / 
            (currentGradientPoint->first - previousGradientPoint->first) * 
            (currentGradientPoint->first - j);
        //std::cout << j << " " << secondSolventConcentration << "<br>";
        S += dV / calculateKd(peptideEnergyProfile, 
                      secondSolventConcentration, 
                      chemBasis, 
                      conditions.columnPoreSize(), 
                      conditions.calibrationParameter(),
                      conditions.temperature()) / volumePore;
        //std::cout << secondSolventConcentration << " " << 
        //1 / (1 + calculateKd (peptideEnergyProfile, 
        //                      chemBasis,
        //                      conditions.secondSolvent(), 
        //                      secondSolventConcentration, 
        //                      conditions.columnPoreSize(), 
        //                      conditions.calibrationParameter(),
        //                      conditions.temperature()) ) << "\n" << "<br>";
    }

    double RT = j * dV / conditions.flowRate() + conditions.delayTime() + 
                volumeLiquidPhase / conditions.flowRate();
    return RT;
}

}

double calculateRT(const std::string &sequence,
    const ChromoConditions &conditions,
    const ChemicalBasis &chemBasis
) {
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;    
    Terminus CTerminus;
    std::vector<double> peptideEnergyProfile;
    
   if (parseSequence(sequence, 
                    chemBasis,
                    &parsedPeptideStructure,
                    &NTerminus,
                    &CTerminus,
                    &peptideEnergyProfile))
    {
        return calculateRT(peptideEnergyProfile,
            conditions,
            chemBasis);
    }
    else {
        return PARSING_ERROR;
    }
}

double calculateKd(const std::string &sequence,
    const double secondSolventConcentration,
    const ChemicalBasis & chemBasis,
    const double columnPoreSize,
    const double calibrationParameter,
    const double temperature
) {
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;    
    Terminus CTerminus;
    std::vector<double> peptideEnergyProfile;
    
    if (parseSequence(sequence, 
                     chemBasis,
                     &parsedPeptideStructure,
                     &NTerminus,
                     &CTerminus,
                     &peptideEnergyProfile))
    {
        return calculateKd(peptideEnergyProfile,
                                  secondSolventConcentration,
                                  chemBasis,
                                  columnPoreSize,
                                  calibrationParameter,
                                  temperature);
    }
    else {
        return PARSING_ERROR;
    }
}                                

double calculateAverageMass(const std::string &sequence,
    const ChemicalBasis &chemBasis
){
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;
    Terminus CTerminus;
    double peptideAverageMass = 0;
    
    if (parseSequence(sequence, chemBasis, &parsedPeptideStructure, 
                      &NTerminus, &CTerminus, NULL)) {
        for (std::vector<Aminoacid>::const_iterator i = 
                 parsedPeptideStructure.begin(); 
             i < parsedPeptideStructure.end(); 
             i++) {
            peptideAverageMass += i -> averageMass();
        }
        peptideAverageMass += NTerminus.averageMass();
        peptideAverageMass += CTerminus.averageMass();
    }
    
    return peptideAverageMass;
}
                                        
double calculateMonoisotopicMass(const std::string &sequence,
    const ChemicalBasis &chemBasis
) {
    std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;
    Terminus CTerminus;
    double monoisotopicMass = 0;
    
    if (parseSequence(sequence, chemBasis, &parsedPeptideStructure, 
                      &NTerminus, &CTerminus, NULL) ) {
        for (std::vector<Aminoacid>::const_iterator i = 
                 parsedPeptideStructure.begin(); 
             i < parsedPeptideStructure.end(); 
             i++) {
            monoisotopicMass += i -> monoisotopicMass();
        }
        monoisotopicMass += NTerminus.monoisotopicMass();
        monoisotopicMass += CTerminus.monoisotopicMass();
    }
    
    return monoisotopicMass;
}


bool calculatePeptideProperties(const std::string &sequence,
         const ChromoConditions &conditions,
         const ChemicalBasis &chemBasis,
         double *RTBioLCCC,
         double *averageMass,
         double *monoisotopicMass
) {
   std::vector<Aminoacid> parsedPeptideStructure;
    Terminus NTerminus;    
    Terminus CTerminus;
    std::vector<double> peptideEnergyProfile;
    
   if (parseSequence(sequence, 
                     chemBasis,
                     &parsedPeptideStructure,
                     &NTerminus,
                     &CTerminus,
                     &peptideEnergyProfile))
    {
        *RTBioLCCC = calculateRT(peptideEnergyProfile, conditions, chemBasis);
        
        *averageMass = 0;
        for (std::vector<Aminoacid>::const_iterator i = 
                 parsedPeptideStructure.begin(); 
             i < parsedPeptideStructure.end(); 
             i++) {
            *averageMass += i -> averageMass();
        }
        *averageMass += NTerminus.averageMass();
        *averageMass += CTerminus.averageMass();
        
        *monoisotopicMass = 0;
        for (std::vector<Aminoacid>::const_iterator i = 
                 parsedPeptideStructure.begin(); 
             i < parsedPeptideStructure.end();
             i++) {
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

}

