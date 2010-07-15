#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
//#include "boost/foreach.hpp"
//#include "boost/function.hpp"
//#include "boost/bind.hpp"
//#include "boost/ref.hpp"
#include "BioLCCC.h"

#define MEMORY_ERROR    -2.0
#define PARSING_ERROR   -3.0
#define PORESIZE_ERROR  -4.0
#define GRADIENT_ERROR  -5.0
#define MODEL_ERROR     -6.0 

#define PI 3.14159265

namespace BioLCCC {

// Auxiliary functions that shouldn't be exposed to user at this
// point.
namespace {
std::vector<double> calculateRT(std::vector<std::string> mixture,
                                ChromoConditions chromatograph,
                                ChemicalBasis chemicalBasis) {
    std::vector<double> times;
    for(unsigned int i = 0; i < mixture.size(); i++) {
        times.push_back(calculateRT(mixture[i], chromatograph, chemicalBasis));
    }
    return times;
}

bool parseSequence(
    const std::string &source, 
    const ChemicalBasis &chemBasis,
    std::vector<Aminoacid> *parsedPeptideStructure,
    Terminus *NTerminus,
    Terminus *CTerminus,
    std::vector<double> *peptideEnergyProfile
) {
    parsedPeptideStructure->clear();
    
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
    //std::pair<std::string,Terminus> NTerminusIterator;
    //BOOST_FOREACH(NTerminusIterator, chemBasis.NTermini()) {
    for (std::map<std::string,Terminus>::const_iterator NTerminusIterator = 
             chemBasis.NTermini().begin();
         NTerminusIterator != chemBasis.NTermini().end();
         NTerminusIterator++) {
        if (strippedSource.find(NTerminusIterator->second.label())== (size_t)0){
            *NTerminus = NTerminusIterator->second;
            NTerminusPosition = NTerminusIterator->second.label().size();
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
        //std::pair<std::string,Terminus> CTerminusIterator;
        //BOOST_FOREACH(CTerminusIterator, chemBasis.CTermini()) {
        for (std::map<std::string,Terminus>::const_iterator CTerminusIterator = 
                 chemBasis.CTermini().begin();
             CTerminusIterator != chemBasis.CTermini().end();
             CTerminusIterator++) {
            if (strippedSource.find(CTerminusIterator->second.label(), 
                CTerminusPosition ) != std::string::npos) {
                *CTerminus = CTerminusIterator->second;
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

    // Then we divide the whole sequence into aminoacids. 
    bool aminoacidFound;
    size_t curPos = 0;
    while (curPos < strippedSource.size()) {
        aminoacidFound = false;
        //std::pair<std::string,Aminoacid> aminoacidIterator;
        //BOOST_FOREACH(aminoacidIterator, chemBasis.aminoacids()) {
        for (std::map<std::string,Aminoacid>::const_iterator aminoacidIterator= 
                 chemBasis.aminoacids().begin();
             aminoacidIterator != chemBasis.aminoacids().end();
             aminoacidIterator++) {
            if (strippedSource.compare(curPos, 
                    aminoacidIterator->second.label().size(), 
                    aminoacidIterator->second.label()) == 0) {
                curPos += aminoacidIterator->second.label().size();
                parsedPeptideStructure -> push_back(aminoacidIterator->second);
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
        //BOOST_FOREACH(Aminoacid currentAminoacid, *parsedPeptideStructure){
        for (std::vector<Aminoacid>::const_iterator currentAminoacid =
                parsedPeptideStructure->begin();
            currentAminoacid != parsedPeptideStructure->end();
            currentAminoacid++) {
            peptideEnergyProfile->push_back(currentAminoacid->bindEnergy());
        }
        
        // Modifing energies of terminal aminoacid residues.
        *(peptideEnergyProfile->begin()) = *(peptideEnergyProfile->begin()) +
            NTerminus->bindEnergy();
        *(--peptideEnergyProfile->end()) = *(--peptideEnergyProfile->end()) + 
            CTerminus->bindEnergy();
    }

    return true;
}

double calculateKdCoilBoltzmann(
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
    Nb = secondSolventConcentration * 1.91 / 
         (secondSolventConcentration * 1.91 + 
          (100.0 - secondSolventConcentration) * 5.56);
    Eab = 0.0 + 1.0 / calibrationParameter * log( 1.0 - Nb + Nb * Q );
   
    // A Boltzmann factor is an exponent of an energy of interaction between
    // an amino acid residue and a solid phase divided by a temperature * 
    // Boltzmann's constant. An energy unit is a Boltzmann's constant * 
    // 293.0 Kelvins. This probability is used later in the transition matrix.
    std::vector<double> boltzmannFactorProfile;
    int residueNumber = 0;
    double energySum = 0.0;
    for (std::vector<double>::const_iterator residueEnergy = 
             peptideEnergyProfile.begin();
         residueEnergy != peptideEnergyProfile.end();
         residueEnergy ++
    ) {
        // More detailed formula for the Nb is:
        // Nb = (DensityB * %B / MB) / (DensityB * %B / MB + DensityA * %A / MA)
        // Where DensityA and DensityB are the corresponding densities and
        // MA and MB are the corresponding molecular weights.
        residueNumber++;
        energySum += calibrationParameter * (*residueEnergy - Eab) *
            293.0 / temperature;
        if ((residueNumber % chemBasis.persistentLength()) == 0) {
            boltzmannFactorProfile.push_back(exp(energySum));
            energySum = 0.0;
        }
    }

    // If the last persistent block is incomplete, we should add it manually.
    if ((residueNumber % chemBasis.persistentLength()) != 0) {
        boltzmannFactorProfile.push_back(exp(energySum));
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
    const int poreSteps = (int) (columnPoreSize / 
        chemBasis.segmentLength() / (double)(chemBasis.persistentLength())) ;
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
        density[i] = 1.0;
    }

    density[poreSteps - 1] = 
        boltzmannFactorProfile[0];

    // Than we construct a basis for the transition matrix. The basis is 
    // a diagonal matrix with 4.0/6.0 on the main diagonal and 1.0/6.0 on 
    // the side diagonals.
    
    // Filling the matrix.
    for (int i = 0; i < poreSteps; i++) {
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

        // Transferring the results from the density vector.
        for (int i = 0; i < poreSteps; i++)  {
            density[i] = densityBuffer[i];
        }
    }

    // Finally, Kd is calculated as a sum of elements of the density vector. 
    // It is normalized to the size of the lattice.
    double Kd=0;
    for (int i=0; i < poreSteps; i++) {
        Kd += density[i];
    }
    Kd = Kd / (double)(poreSteps);

    // Cleaning memory.
    try {
        delete[] density;
        delete[] densityBuffer;
        delete[] transitionMatrix;
    }
    catch (...) {
        return MEMORY_ERROR;
    }

    return Kd;
}                                            

double calculateKdCoilBoltzmannDoubleLayer(
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
    // the features of the following calculation it should be more than 2.
    const int poreSteps = (int) (columnPoreSize / chemBasis.segmentLength());
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
    density[1] = boltzmannFactorProfile[0];

    for (int i = 2; i < poreSteps - 2 ; i++) {
        density[i] = 1;
    }

    density[poreSteps - 2] = 
        boltzmannFactorProfile[0];
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
        transitionMatrix[poreSteps] = 1.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[poreSteps + 1] = 4.0 / 6.0 * (*residueBoltzmannFactor);
        transitionMatrix[poreSteps + 2] = 1.0 / 6.0 * (*residueBoltzmannFactor);

        transitionMatrix[poreSteps*(poreSteps - 1) - 1] = 1.0 / 6.0 * 
            (*residueBoltzmannFactor);
        transitionMatrix[poreSteps*(poreSteps - 1) - 2] = 4.0 / 6.0 * 
            (*residueBoltzmannFactor);
        transitionMatrix[poreSteps*(poreSteps - 1) - 3] = 1.0 / 6.0 * 
            (*residueBoltzmannFactor);
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
    for (int k=0; k < poreSteps; k++) {
        Kd += density[k];
    }
    Kd = Kd / poreSteps;

    // Cleaning memory.
    try {
        delete[] density;
        delete[] densityBuffer;
        delete[] transitionMatrix;
    }
    catch (...) {
        return MEMORY_ERROR;
    }

    return Kd;
}                                            

double calculateKdCoilSnyder(
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
    double Nb = 0;
   
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
        boltzmannFactorProfile.push_back(
            exp(calibrationParameter * 
                (*residueEnergy - chemBasis.secondSolventBindEnergy() * Nb) *
                293.0 / temperature));
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
    const int poreSteps = (int)(columnPoreSize / chemBasis.segmentLength());
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
        delete[] density;
        delete[] densityBuffer;
        delete[] transitionMatrix;
    }
    catch (...) {
        return MEMORY_ERROR;
    }

    return Kd;
}                                            

double rodAdsorbtionEnergy(const std::vector<double> & peptideEnergyProfile,
                           unsigned int n,
                           bool reversed = false
) {
    double init = 0;
    if (n > peptideEnergyProfile.size()){
        // TODO: make an exception.
        return -65536;
    }
    if (reversed) {
        return std::accumulate(peptideEnergyProfile.end()-n,
                               peptideEnergyProfile.end(),
                               init);
    }
    else {
        return std::accumulate(peptideEnergyProfile.begin(),
                               peptideEnergyProfile.begin()+n,
                               init);
    }
}

double partitionFunctionRodFreeSlit(double rodLength,
                                    double slitWidth
) {
    // the equation works only if the slit is wider than the rod
    if (rodLength <= slitWidth) {
    // full volume without exclusion caused by walls
        return (4 * PI * slitWidth * rodLength * rodLength -
    // minus volume excluded by walls
                2 * PI * rodLength * rodLength * rodLength);
    }
    else {
    // This case is considered in the paper as Zc.
        return (2 * PI * slitWidth * slitWidth * rodLength);
    }
        
}

double partitionFunctionRodFreeSlitA(double rodLength,
                                     double slitWidth
) {
        return (4 * PI * slitWidth * rodLength * rodLength -
                2 * PI * rodLength * rodLength * rodLength);
}

double partitionFunctionRodFreeSlitC(double rodLength,
                                     double slitWidth
) {
        return (2 * PI * slitWidth * slitWidth * rodLength);
}

double partitionFunctionRodSubmergedIntoLayer(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed = false
) {
    double rodLength = segmentLength * (double)(rodEnergyProfile.size() - 1);
    double partitionFunction = 0.0;
    for (unsigned int n = 1; n < rodEnergyProfile.size(); n++) {
        if (layerWidth >= segmentLength * double(n)) {
            partitionFunction += 2.0 * PI * rodLength * rodLength *
                exp(rodAdsorbtionEnergy(rodEnergyProfile, n, reversed)) *
                segmentLength / 2.0;
        }
        else if ((segmentLength * (double)(n-1) < layerWidth) && 
                 (layerWidth < segmentLength * double(n))) {
            partitionFunction += 2.0 * PI * rodLength * rodLength *
                exp(rodAdsorbtionEnergy(rodEnergyProfile, n, reversed)) *
                (layerWidth
                 - segmentLength * (double)(n-1) / 2.0
                 - layerWidth * layerWidth / 2.0 / (double)n / segmentLength);
        }
        else {
            partitionFunction += 2.0 * PI * rodLength * rodLength *
                exp(rodAdsorbtionEnergy(rodEnergyProfile, n, reversed)) *
                layerWidth * layerWidth / 2.0 / double(n) / double(n-1) /
                segmentLength;
        }
    }
    return partitionFunction;
}

double partitionFunctionRodFreeVolume(double rodLength,
                                      double slitWidth
) {
    return (4 * PI * slitWidth * rodLength * rodLength);
}

double calculateKdRod(
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
    std::vector<double> effectiveEnergyProfile;

    //BOOST_FOREACH(double residueEnergy, peptideEnergyProfile) {
    for (std::vector<double>::const_iterator residueEnergy=
            peptideEnergyProfile.begin();
        residueEnergy != peptideEnergyProfile.end();
        residueEnergy++) {
        Nb = secondSolventConcentration * 1.91 / 
             (secondSolventConcentration * 1.91 + 
              (100.0 - secondSolventConcentration) * 5.56);
        Eab = 0.0 + 1.0 / calibrationParameter * log( 1.0 - Nb + Nb * Q );
        effectiveEnergyProfile.push_back(calibrationParameter * 
            (*residueEnergy - Eab) * 293.0 / temperature);
    }
   
    //std::cout << "Calculating Kd:\n";
    //std::cout << 
    //    partitionFunctionRodFreeSlit(
    //    chemBasis.segmentLength() * (effectiveEnergyProfile.size() - 1),
    //    columnPoreSize)<< "\n";

    //std::cout << 
    //    2.0 * partitionFunctionRodFreeSlit(
    //            chemBasis.segmentLength() * (effectiveEnergyProfile.size() - 1),
    //            chemBasis.adsorbtionLayerWidth()) 
    //    * exp(rodAdsorbtionEnergy(
    //        effectiveEnergyProfile,
    //        effectiveEnergyProfile.size()))<< "\n";

    //std::cout << 
    //2.0 * partitionFunctionRodSubmergedIntoLayer(
    //    chemBasis.segmentLength(),
    //    columnPoreSize,
    //    chemBasis.adsorbtionLayerWidth(),
    //    effectiveEnergyProfile,
    //    false) << "\n";

    //std::cout << 2.0 * partitionFunctionRodSubmergedIntoLayer(
    //    chemBasis.segmentLength(),
    //    columnPoreSize,
    //    chemBasis.adsorbtionLayerWidth(),
    //    effectiveEnergyProfile,
    //    true) << "\n";

    //std::cout << 
    //    partitionFunctionRodFreeVolume(
    //    (effectiveEnergyProfile.size() - 1) * chemBasis.segmentLength(),
    //    columnPoreSize) << "\n"; 

    double Kd = 
    (partitionFunctionRodFreeSlit(
        chemBasis.segmentLength() * (effectiveEnergyProfile.size() - 1),
        columnPoreSize - 2.0 * chemBasis.adsorbtionLayerWidth()) 
    + 2.0 * partitionFunctionRodFreeSlit(
                chemBasis.segmentLength() * (effectiveEnergyProfile.size() - 1),
                chemBasis.adsorbtionLayerWidth()) 
        * exp(rodAdsorbtionEnergy(
            effectiveEnergyProfile,
            effectiveEnergyProfile.size()))

    + 2.0 * partitionFunctionRodSubmergedIntoLayer(
        chemBasis.segmentLength(),
        columnPoreSize,
        chemBasis.adsorbtionLayerWidth(),
        effectiveEnergyProfile,
        false) 
    + 2.0 * partitionFunctionRodSubmergedIntoLayer(
        chemBasis.segmentLength(),
        columnPoreSize,
        chemBasis.adsorbtionLayerWidth(),
        effectiveEnergyProfile,
        true)) /
    partitionFunctionRodFreeVolume(
        (effectiveEnergyProfile.size() - 1) * chemBasis.segmentLength(),
        columnPoreSize); 

    return Kd;
}                                            

double calculateKd(const std::vector<double> &peptideEnergyProfile,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double calibrationParameter,
    const double temperature
) {
    // Choosing the appropriate model.
    if (chemBasis.model()=="CoilBoltzmann") {
        if (peptideEnergyProfile.size() % chemBasis.persistentLength() == 0) {
            return calculateKdCoilBoltzmann(peptideEnergyProfile,
                secondSolventConcentration, chemBasis, columnPoreSize,
                calibrationParameter, temperature);
        }
        else {
            std::vector<double> revPeptideEnergyProfile = peptideEnergyProfile;
            std::reverse(revPeptideEnergyProfile.begin(),
                         revPeptideEnergyProfile.end());
            return (calculateKdCoilBoltzmann(peptideEnergyProfile,
                secondSolventConcentration, chemBasis, columnPoreSize,
                calibrationParameter, temperature) +
                calculateKdCoilBoltzmann(revPeptideEnergyProfile,
                secondSolventConcentration, chemBasis, columnPoreSize,
                calibrationParameter, temperature)) / 2.0 ;
        }
    }
    else if (chemBasis.model()=="CoilBoltzmannDoubleLayer") {
        return calculateKdCoilBoltzmannDoubleLayer(peptideEnergyProfile,
            secondSolventConcentration,
            chemBasis,
            columnPoreSize,
            calibrationParameter,
            temperature);
    }
    else if (chemBasis.model() == "CoilSnyder") {
        return calculateKdCoilSnyder(peptideEnergyProfile,
            secondSolventConcentration,
            chemBasis,
            columnPoreSize,
            calibrationParameter,
            temperature);
    }
    else if (chemBasis.model() == "RodBoltzmann") {
        return calculateKdRod(peptideEnergyProfile,
            secondSolventConcentration,
            chemBasis,
            columnPoreSize,
            calibrationParameter,
            temperature);
    }
    else {
        return MODEL_ERROR;
    }
}                                

double calculateRT(const std::vector<double> &peptideEnergyProfile,
    const ChromoConditions &conditions,
    const ChemicalBasis &chemBasis,
    const bool continueGradient
) {
    // Calculating column volumes.
    double volumeLiquidPhase = conditions.columnDiameter() * 
        conditions.columnDiameter() * 3.1415 * conditions.columnLength() / 4.0 /
        1000.0 * (conditions.columnPorosity()-conditions.columnVpToVtot());
    double volumePore = conditions.columnDiameter() * 
        conditions.columnDiameter() * 3.1415 * conditions.columnLength() / 4.0 /
        1000.0 * conditions.columnVpToVtot();    
    
    // Recalculating dV. By default dV is calculated as flow rate divided by 20.
    double dV;
    if (conditions.dV()!= 0.0) {
        dV = conditions.dV();
    }
    else {
        dV = conditions.flowRate() / 20.0;
    }

    // Because of the features of the BioLCCC model, size of the pore should be 
    // more than 20A.
    if (conditions.columnPoreSize() <= 15.0) {
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
    
    // We've dropped the const_iterator version, since it didn't work 
    // on WinGW32.
    //BOOST_FOREACH(GradientPoint currentGradientPoint, conditions.gradient()){
    //for (Gradient::const_iterator currentGradientPoint =
    //        conditions.gradient().begin();
    //    currentGradientPoint < conditions.gradient().end();
    //    ++currentGradientPoint) {
    for (Gradient::size_type i = 0;
        i != conditions.gradient().size();
        i++) {
        convertedGradient.push_back(
                std::pair<int,double>(
                    int(floor(conditions.gradient()[i].time() * 
                         conditions.flowRate() / dV)),
                    (100.0 - conditions.gradient()[i].concentrationB())/100.0 *
                        conditions.secondSolventConcentrationA() + 
                        conditions.gradient()[i].concentrationB() / 100.0 *
                        conditions.secondSolventConcentrationB()));
    }
    
    // The part of a column passed by molecules. When it exceeds 1.0,
    // molecule elute from the column.
    double S = 0.0;
    // The current iteration number. 
    int j = 0;      
    double secondSolventConcentration = 0.0;


    std::vector<std::pair<int, double> >::const_iterator currentGradientPoint=
        convertedGradient.begin();
    std::vector<std::pair<int, double> >::const_iterator previousGradientPoint=
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
                        calculateKd(
                            peptideEnergyProfile, 
                            currentGradientPoint->second, 
                            chemBasis, 
                            conditions.columnPoreSize(), 
                            conditions.calibrationParameter(),
                            conditions.temperature()) *
                        volumePore < (currentGradientPoint->first - j + 1)) ||
                        
                        (currentGradientPoint == --convertedGradient.end()))
                    {
                        j += (int)ceil((1.0 - S) / dV *
                            calculateKd(
                                peptideEnergyProfile, 
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
                        S += dV / calculateKd(
                                      peptideEnergyProfile,
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
            else if (!continueGradient) {
                break;
            }
        }
        secondSolventConcentration = currentGradientPoint->second - 
            (currentGradientPoint->second - previousGradientPoint->second) / 
            (currentGradientPoint->first - previousGradientPoint->first) * 
            (currentGradientPoint->first - (double) j);
        S += dV / calculateKd(
                      peptideEnergyProfile, 
                      secondSolventConcentration, 
                      chemBasis, 
                      conditions.columnPoreSize(), 
                      conditions.calibrationParameter(),
                      conditions.temperature()) / volumePore;
    }

    double RT = j * dV / conditions.flowRate() + conditions.delayTime() + 
                volumeLiquidPhase / conditions.flowRate();
    return RT;
}

}

double calculateRT(const std::string &sequence,
    const ChromoConditions &conditions,
    const ChemicalBasis &chemBasis,
    const bool continueGradient
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
            chemBasis,
            continueGradient);
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
        *RTBioLCCC = calculateRT(peptideEnergyProfile, conditions, chemBasis,
                                 true);
        
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

//ChemicalBasis calibrateBioLCCC(std::vector<std::string> calibrationMixture,
//                               std::vector<double> retentionTimes,
//                               ChromoConditions chromatograph,
//                               ChemicalBasis initialChemicalBasis,
//                               std::vector<std::string> energiesToCalibrate
//) {
//
//    // Generate the list of setter functions out of 'energiesToCalibrate'
//    ChemicalBasis newChemicalBasis = initialChemicalBasis;
//    std::vector< boost::function<void(double)> > setters;
//    for(unsigned int i = 0; i < energiesToCalibrate.size(); i++) {
//        if (initialChemicalBasis.aminoacids().find(energiesToCalibrate[i]) !=
//        initialChemicalBasis.aminoacids().end()) {
//            setters.push_back(boost::bind(
//                &ChemicalBasis::setAminoacidBindEnergy,
//                &newChemicalBasis,
//                energiesToCalibrate[i],
//                _1));
////            std::cout << "AA found\n";
//        }
//        else if (initialChemicalBasis.NTermini().find(energiesToCalibrate[i]) !=
//                 initialChemicalBasis.NTermini().end()) {
//                     setters.push_back(boost::bind(
//                         &ChemicalBasis::setNTerminusBindEnergy,
//                         &newChemicalBasis,
//                         energiesToCalibrate[i],
//                         _1));
////                     std::cout << "NT found\n";
//                 }
//        else if (initialChemicalBasis.CTermini().find(energiesToCalibrate[i]) !=
//                 initialChemicalBasis.CTermini().end()) {
//                     setters.push_back(boost::bind(
//                         &ChemicalBasis::setCTerminusBindEnergy,
//                         &newChemicalBasis,
//                         energiesToCalibrate[i],
//                         _1));
////                     std::cout << "CT found\n";
//                 }
//        else if (energiesToCalibrate[i] == "ACN") {
//            // This is a walkaround for the strange type-punning GCC warning.
//            boost::function<void (ChemicalBasis*, double)> secondSolventSetter;
//            secondSolventSetter = &ChemicalBasis::setSecondSolventBindEnergy;
//            setters.push_back(boost::bind(secondSolventSetter,
//                                          &newChemicalBasis,
//                                          _1));
////            std::cout << "ACN found\n";
//        }
//
//    }
//
//    std::vector<double> (*calculateRTforList)
//        (std::vector<std::string>,ChromoConditions,ChemicalBasis) =
//            &calculateRT;
//    boost::function<double(void)> penaltyFunction = boost::bind(
//        boost::bind(vectorNorm<double>, 
//                    boost::bind(VectorDiff<double>, _1, _2)),
//        boost::bind(calculateRTforList, calibrationMixture, 
//                    chromatograph, boost::ref(newChemicalBasis)),
//        retentionTimes);
//
//    std::vector<double> lowerBounds, upperBounds, steps, gradientSteps;
//    for(unsigned int i = 0; i < setters.size(); i++) {
//        lowerBounds.push_back(0.0);
//        upperBounds.push_back(3.0);
//        steps.push_back(0.25);
//        gradientSteps.push_back(0.001);
//    }
//    //std::vector<double> lowerBounds, upperBounds, steps, gradientSteps;
//    //lowerBounds.push_back(0.575-0.01);
//    //lowerBounds.push_back(0.83-0.01);
//    //lowerBounds.push_back(2.20-0.01);
//    //upperBounds.push_back(0.575+0.01);
//    //upperBounds.push_back(0.83+0.01);
//    //upperBounds.push_back(2.20+0.01);
//    //steps.push_back(0.0025);
//    //steps.push_back(0.0025);
//    //steps.push_back(0.0025);
//    //gradientSteps.push_back(0.001);
//    //gradientSteps.push_back(0.001);
//    //gradientSteps.push_back(0.001);
//
//    std::cout << "Initial penalty function value: " 
//        << penaltyFunction() << "\n";
//    std::vector<double> calibratedEnergiesBruteForce = 
//        findMinimumBruteForce(penaltyFunction, setters, lowerBounds,
//                              upperBounds, steps);
//    std::cout << "BruteForce approximation:\n";
//    for(unsigned int i = 0; i < setters.size(); i++) {
//        std::cout << energiesToCalibrate[i] << ": " 
//            << calibratedEnergiesBruteForce[i] << "\n";
//    }
//
//    for(unsigned int i = 0; i < setters.size(); i++) {
//        setters[i](calibratedEnergiesBruteForce[i]);
//    }
//    std::cout << "Penalty function value after calibration: " 
//        << penaltyFunction() << "\n";
//
//    std::vector<double> calibratedEnergiesGradientDescent = 
//        findMinimumGradientDescent(penaltyFunction, setters, 
//            calibratedEnergiesBruteForce, gradientSteps, 0.00001);
//
//    std::cout << "GradientDescent approximation:\n";
//    for(unsigned int i = 0; i < setters.size(); i++) {
//        std::cout << energiesToCalibrate[i] << ": " 
//            << calibratedEnergiesGradientDescent[i] << "\n";
//    }
//
//    for(unsigned int i = 0; i < setters.size(); i++) {
//        setters[i](calibratedEnergiesGradientDescent[i]);
//    }
//    return newChemicalBasis;
//}

}
