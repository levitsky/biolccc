#include <cmath>
#include <numeric>
#include "rod_model.h"
#include "parsing.h"

#define PI 3.14159265

namespace BioLCCC
{

double rodAdsorptionEnergy(const std::vector<double> & rodEnergyProfile,
                           unsigned int n,
                           bool reversed) throw(BioLCCCException)
{
    if ((n < 0) || (n > rodEnergyProfile.size()))
    {
        throw BioLCCCException("Index is out of range.");
    }

    double init = 0;
    if (reversed)
    {
        return std::accumulate(rodEnergyProfile.end()-n,
                               rodEnergyProfile.end(),
                               init);
    }
    else
    {
        return std::accumulate(rodEnergyProfile.begin(),
                               rodEnergyProfile.begin()+n,
                               init);
    }
}

double partitionFunctionRodFreeSlit(double rodLength,
                                    double slitWidth)
{
    // the equation works only if the slit is wider than the rod
    if (rodLength <= slitWidth)
    {
        // full volume without exclusion caused by walls
        return (4 * PI * slitWidth * rodLength * rodLength -
                // minus volume excluded by walls
                2 * PI * rodLength * rodLength * rodLength);
    }
    else
    {
        // This case is considered in the paper as Zc.
        return (2 * PI * slitWidth * slitWidth * rodLength);
    }

}

//double partitionFunctionRodSubmergedIntoLayerGeneral(
//    double segmentLength,
//    double slitWidth,
//    double layerWidth,
//    const std::vector<double> & rodEnergyProfile,
//    bool reversed = false) throw(BioLCCCException)
//{
//    double partitionFunction = 0.0;
//    const double N = rodEnergyProfile.size();
//    const double rodLength = segmentLength * (N - 1);
//    for (n1 = 1; n1 < N; n1++)
//    {
//        for (n2 = n1 + 1; n2 <= N; n2++)
//        {
//            partitionFunction +=
//                
//                * exp(rodAdsorptionEnergy(rodEnergyProfile, n, reversed));
//        }
//    }
//}

double partitionFunctionRodSubmergedIntoLayer(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed) throw(BioLCCCException)
{
    double rodLength = segmentLength * (double)(rodEnergyProfile.size() - 1);
    double partitionFunction = 0.0;
    for (unsigned int n = 1; n < rodEnergyProfile.size(); n++)
    {
        if (layerWidth >= segmentLength * double(n))
        {
            partitionFunction += 2.0 * PI * rodLength * rodLength *
                                 exp(rodAdsorptionEnergy(
                                    rodEnergyProfile, n, reversed)) *
                                 segmentLength / 2.0;
        }
        else if ((segmentLength * (double)(n-1) < layerWidth) &&
                 (layerWidth < segmentLength * double(n)))
        {
            partitionFunction += 2.0 * PI * rodLength * rodLength *
                                 exp(rodAdsorptionEnergy(
                                    rodEnergyProfile, n, reversed)) *
                                 (layerWidth
                                  - segmentLength * (double)(n-1) / 2.0
                                  - layerWidth * layerWidth / 2.0 / 
                                      (double)n / segmentLength);
        }
        else
        {
            partitionFunction += 2.0 * PI * rodLength * rodLength *
                                 exp(rodAdsorptionEnergy(
                                    rodEnergyProfile, n, reversed)) *
                                 layerWidth * layerWidth / 2.0 / double(n)
                                 / double(n-1) / segmentLength;
        }
    }
    return partitionFunction;
}

double partitionFunctionRodFreeVolume(double rodLength,
                                      double slitWidth)
{
    return (4 * PI * slitWidth * rodLength * rodLength);
}

double calculateKdRod(
    const std::vector<ChemicalGroup> &parsedSequence,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double columnRelativeStrength,
    const double temperature) throw(BioLCCCException)
{
    if (parsedSequence.size() == 0)
    {
        return 0.0;
    }

    std::vector<double> segmentEnergyProfile = 
        calculateSegmentEnergyProfile(
            calculateMonomerEnergyProfile(
                parsedSequence,
                chemBasis,
                secondSolventConcentration,
                columnRelativeStrength,
                temperature),
            chemBasis.monomerLength(),
            chemBasis.kuhnLength());

    double rodLength = chemBasis.kuhnLength() *
        (segmentEnergyProfile.size() - 1);

    double Kd =
        (
            partitionFunctionRodFreeSlit(
                rodLength,
                columnPoreSize - 2.0 * chemBasis.adsorptionLayerWidth())

            + 2.0 * partitionFunctionRodFreeSlit(
                rodLength,
                chemBasis.adsorptionLayerWidth())
            * exp(rodAdsorptionEnergy(
                segmentEnergyProfile,
                segmentEnergyProfile.size()))

            + 2.0 * partitionFunctionRodSubmergedIntoLayer(
                chemBasis.kuhnLength(),
                columnPoreSize,
                chemBasis.adsorptionLayerWidth(),
                segmentEnergyProfile,
                false)

            + 2.0 * partitionFunctionRodSubmergedIntoLayer(
                chemBasis.kuhnLength(),
                columnPoreSize,
                chemBasis.adsorptionLayerWidth(),
                segmentEnergyProfile,
                true)
        ) / partitionFunctionRodFreeVolume(
            rodLength,
            columnPoreSize);

    return Kd;
}

}
