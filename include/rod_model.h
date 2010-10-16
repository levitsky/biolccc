#ifndef ROD_MODEL_H
#define ROD_MODEL_H

#include "biolcccexception.h"
#include "chemicalbasis.h"

namespace BioLCCC
{

//! Calculates the adsorption energy of first n segments of the rod.
double rodAdsorptionEnergy(const std::vector<double> & rodEnergyProfile,
                           unsigned int n,
                           bool reversed = false) throw(BioLCCCException);

//! Calculates the partition function of the rod in a slit.
double partitionFunctionRodFreeSlit(double rodLength,
                                    double slitWidth);

//! Calculates Z of the rod partially submerged into the adsorbing layer.
double partitionFunctionRodSubmergedIntoLayer(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed = false) throw(BioLCCCException);

//! Calculates the partition function of the rod in a slit ignoring the walls.
double partitionFunctionRodFreeVolume(double rodLength,
                                      double slitWidth);

//! Calculates the coefficient of distribution of a polymer using the rod model.
double calculateKdRod(
    const std::vector<ChemicalGroup> &parsedSequence,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double columnRelativeStrength,
    const double temperature) throw(BioLCCCException);

}

#endif
