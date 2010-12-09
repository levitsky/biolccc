#ifndef ROD_MODEL_H
#define ROD_MODEL_H

#include "biolcccexception.h"
#include "chemicalbasis.h"

namespace BioLCCC
{

double partitionFunctionRodPartiallySubmergedTermSpecial(
    double segmentLength, double slitWidth, double layerWidth,
    int N, int n1);

double partitionFunctionRodPartiallySubmergedTermGeneral(
    double segmentLength, double slitWidth, double layerWidth,
    int N, int n1, int n2);

//! Calculates the adsorption energy of first n segments of the rod.
double rodAdsorptionEnergy(const std::vector<double> & rodEnergyProfile,
                           int n1, int n2) throw(BioLCCCException);

//! Calculates the partition function of the rod in a slit.
double partitionFunctionRodFreeSlit(double rodLength,
                                    double slitWidth);

//! Calculates Z of the rod partially submerged into the adsorbing layer.
double partitionFunctionRodPartiallySubmergedGeneral(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed = false) throw(BioLCCCException);

//! Calculates Z of the rod partially submerged into the adsorbing layer.
double partitionFunctionRodPartiallySubmergedSpecial(
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
    const double temperature,
    const bool specialRodModel
    ) throw(BioLCCCException);
}

#endif
