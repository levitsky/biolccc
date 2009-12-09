#ifndef PEPTIDEMETHODS_H
#define PEPTIDEMETHODS_H

#include <math.h>
#include "chemicalbasis.h"
#include "chromoconditions.h"
#include "auxiliary.hpp"

//! This is the main header for the whole BioLCCC library project. 

/*!
    <b> Known issues </b>
    - Due to the fundamental ambiguity you couldn't parse properly the 
    sequence of amidated histidine (H-NH2).
*/

namespace BioLCCC {

static const ChromoConditions standardChromoConditions = ChromoConditions();
static const ChemicalBasis standardChemicalBasis = ChemicalBasis();

/*!
    Calculates the retention time of a peptide with given sequence 
    using given table of peptide chemicals and chromatographic 
    conditions.
*/
double calculateRT(const std::string &sequence,
    const ChromoConditions & conditions = standardChromoConditions,
    const ChemicalBasis & chemBasis = standardChemicalBasis,
    const bool continueGradient = true);

/*!
    Calculates the average (molar) mass of a peptide with given 
    sequence using given table of peptide chemicals.
*/
double calculateAverageMass(const std::string &sequence,
    const ChemicalBasis &chemBasis = standardChemicalBasis);
                        
/*!
    Calculates the monoisotopic mass of a peptide with given sequence 
    using given table of peptide chemicals.
*/
double calculateMonoisotopicMass(const std::string &sequence,
    const ChemicalBasis &chemBasis = standardChemicalBasis);

/*!
    Calculates the distribution coefficient (Kd) of a peptide with
    given sequence with the BioLCCC model using given table of peptide 
    chemicals, the name of the second solvent, its concentration, 
    the size of adsorbent's pores, calibration parameter and 
    temperature.
*/
double calculateKd(const std::string &sequence,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis = standardChemicalBasis,
    const double columnPoreSize = 100.0,
    const double calibrationParameter = 1.0,
    const double temperature = 293.0);

/*!
    Created as a transient solution for fast calculation of RTBioLCCC 
    and masses.
*/
bool calculatePeptideProperties(const std::string &sequence,
    const ChromoConditions &conditions,
    const ChemicalBasis &chemBasis,
    double *RTBioLCCC,
    double *averageMass,
    double *monoisotopicMass);

ChemicalBasis calibrateBioLCCC(
        std::vector<std::string> calibrationMixture,
        std::vector<double> retentionTimes,
        ChromoConditions chromatograph,
        ChemicalBasis initialChemicalBasis,
        std::vector<std::string> energiesToCalibrate);

}
#endif
