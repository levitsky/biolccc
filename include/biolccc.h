#ifndef BIOLCCC_H
#define BIOLCCC_H

#include <math.h>
#include "biolcccexception.h"
#include "chemicalbasis.h"
#include "chromoconditions.h"

//! This is the main header for the whole BioLCCC library project.

/*!
    <b> Known issues </b>
    - Due to the fundamental ambiguity you couldn't parse properly the
    sequence of amidated histidine (H-NH2).
*/

namespace BioLCCC
{

class ParsingException : public BioLCCCException
{
public:
    ParsingException(std::string message);
};

static const ChromoConditions standardChromoConditions = ChromoConditions();
static const ChemicalBasis standardChemicalBasis = ChemicalBasis();

/*!
    Parses a given string representation of a peptide using given chemical
    basis. Returns parsed peptide structure, terminal groups and energy profile
    of a peptide.
*/
void parseSequence(
    const std::string &source,
    const ChemicalBasis &chemBasis,
    std::vector<ChemicalGroup> *parsedPeptideStructure,
    ChemicalGroup *NTerminus,
    ChemicalGroup *CTerminus,
    std::vector<double> *peptideEnergyProfile);

/*!
    Calculates the retention time of a peptide with given sequence
    using given table of peptide chemicals and chromatographic
    conditions.
*/
double calculateRT(const std::string &sequence,
                   const ChromoConditions & conditions = 
                       standardChromoConditions,
                   const ChemicalBasis & chemBasis = standardChemicalBasis,
                   const bool continueGradient = true);

/*!
    Calculates the average (molar) mass of a peptide with given
    sequence using given table of peptide chemicals.
*/
double calculateAverageMass(const std::string &sequence,
                            const ChemicalBasis &chemBasis = 
                                standardChemicalBasis);

/*!
    Calculates the monoisotopic mass of a peptide with given sequence
    using given table of peptide chemicals.
*/
double calculateMonoisotopicMass(const std::string &sequence,
                                 const ChemicalBasis &chemBasis = 
                                     standardChemicalBasis);

/*!
    Calculates the distribution coefficient (Kd) of a peptide with
    given sequence with the STANDARD BioLCCC model (for a coil using the
    Boltzmann equation) with given table of peptide chemicals, the name of
    the second solvent, its concentration, the size of adsorbent's pores,
    calibration parameter and temperature.
*/
double calculateKd(const std::string &sequence,
                   const double secondSolventConcentration,
                   const ChemicalBasis &chemBasis = standardChemicalBasis,
                   const double columnPoreSize = 100.0,
                   const double calibrationParameter = 1.0,
                   const double temperature = 293.0);

/*!
    Created as a transient solution for the fast calculation of RTBioLCCC
    and masses.
*/
void calculatePeptideProperties(const std::string &sequence,
                                const ChromoConditions &conditions,
                                const ChemicalBasis &chemBasis,
                                double *RTBioLCCC,
                                double *averageMass,
                                double *monoisotopicMass);

}
#endif
