#ifndef BIOLCCC_H
#define BIOLCCC_H

#include <math.h>
#include "biolcccexception.h"
#include "chemicalbasis.h"
#include "chromoconditions.h"

namespace BioLCCC
{

//! This exception is raised when parsing process cannot be completed.
class ParsingException : public BioLCCCException
{
public:
    //! Constructs an instance ParsingException with the given \a message.
    ParsingException(std::string message);
};


static const ChromoConditions standardChromoConditions = ChromoConditions();
static const ChemicalBasis standardChemicalBasis = ChemicalBasis();

//! Parses the given peptide sequence.
/*!
    Parses the given peptide sequence \a source using \a chemBasis. Writes
    the parsed peptide structure into \a parsedPeptideStructure, terminal groups
    into \a NTerminus and \a CTerminus. Writes the energy profile of a peptide
    into \a peptideEnergyProfile.

    Throws ParsingException if the peptide is not parsable.
*/
void parseSequence(
    const std::string &source,
    const ChemicalBasis &chemBasis,
    std::vector<ChemicalGroup> *parsedPeptideStructure,
    ChemicalGroup *NTerminus,
    ChemicalGroup *CTerminus,
    std::vector<double> *peptideEnergyProfile);

//! Calculates the retention time of a peptide.
/*!
    Calculates the retention time of a peptide with given \a sequence
    using the given description of chromatographic conditions \a conditions and
    set of physicochemical constants \a chemBasis.
    
    If \a continueGradient is true, than the last section of a gradient is
    prolonged.
*/
double calculateRT(const std::string &sequence,
                   const ChromoConditions & conditions = 
                       standardChromoConditions,
                   const ChemicalBasis & chemBasis = standardChemicalBasis,
                   const bool continueGradient = true);

//! Calculates the average (molar) mass of a peptide.
/*!
    Calculates the average (molar) mass of a peptide with given
    \a sequence using the given set of physicochemical constants \a chemBasis.
*/
double calculateAverageMass(const std::string &sequence,
                            const ChemicalBasis &chemBasis = 
                                standardChemicalBasis);

//! Calculates the monoisotopic mass of a peptide.
/*!
    Calculates the monoisotopic mass of a peptide with given
    \a sequence using the given set of physicochemical constants \a chemBasis.
*/
double calculateMonoisotopicMass(const std::string &sequence,
                                 const ChemicalBasis &chemBasis = 
                                     standardChemicalBasis);

//! Calculates the coefficient of distribution Kd for the given peptide.
/*!
    Calculates the coefficient of distribution Kd (i.e. the ratio of
    concentrations of a peptide in the pores and in the interstitial volume).
    
    \param sequence The sequence of a peptide.
    \param secondSolventConcentration The concentration of the second solvent in
    the liquid phase
    \param chemBasis The set of the physicochemical constants.
    \param columnPoreSize The size of adsorbent pores.
    \param columnRelativeStrength The relative strength of adsorption.
    \param temperature Temperature of the column.
*/
double calculateKd(const std::string &sequence,
                   const double secondSolventConcentration,
                   const ChemicalBasis &chemBasis = standardChemicalBasis,
                   const double columnPoreSize = 100.0,
                   const double columnRelativeStrength = 1.0,
                   const double temperature = 293.0);

///*!
//    Created as a transient solution for the fast calculation of RTBioLCCC
//    and masses.
//*/
//void calculatePeptideProperties(const std::string &sequence,
//                                const ChromoConditions &conditions,
//                                const ChemicalBasis &chemBasis,
//                                double *RTBioLCCC,
//                                double *averageMass,
//                                double *monoisotopicMass);
//
}
#endif
