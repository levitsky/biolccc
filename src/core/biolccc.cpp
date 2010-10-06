#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "biolccc.h"

#define PI 3.14159265

namespace BioLCCC
{

ParsingException::ParsingException(std::string message):
        BioLCCCException(message) 
{
};

// Auxiliary functions that shouldn't be exposed to a user.
namespace
{

/*!
    This function calculates the effective energy profile of a polymer chain.
    Each element in this profile contains the adsorption energy of a single
    monomer.

    E_{effective} = alpha * ( E_{residue} - E_{ab} ) * 293.0 / temperature,
    and Eab is the effective energy of binary solvent binding to the
    stationary phase.
    Eab = Ea + 1 / alpha * ln( 1 + Nb + Nb * exp( alpha * ( Ea - Eb ) ) )
*/
std::vector<double> calculateMonomerEnergyProfile(
    const std::vector<ChemicalGroup> &parsedSequence,
    const ChemicalBasis & chemBasis,
    const double secondSolventConcentration,
    const double columnRelativeStrength, 
    const double temperature) throw (BioLCCCException)
{
    if (parsedSequence.size() < 3)
    {
        throw BioLCCCException(
            "The parsed sequence contains too little chemical groups.");
    }

    if (columnRelativeStrength == 0.0)
    {
        return std::vector<double> (parsedSequence.size()-2, 0.0);
    }

    // Due to the preliminary scaling the binding energy of water equals zero.
    double Q = exp(columnRelativeStrength * 
                   (chemBasis.secondSolventBindEnergy() - 0.0) *
                   293.0 / temperature);
    // Nb = (DensityB * %B / MB) / 
    //      (DensityB * %B / MB + DensityA * %A / MA)
    // where DensityA and DensityB are the corresponding densities and
    // MA and MB are the corresponding molecular weights.
    double Nb =
        secondSolventConcentration * chemBasis.secondSolventDensity()
        / chemBasis.secondSolventAverageMass() 
        / ( secondSolventConcentration * chemBasis.secondSolventDensity() 
                / chemBasis.secondSolventAverageMass()
            + (100.0 - secondSolventConcentration) 
              * chemBasis.firstSolventDensity()
              / chemBasis.firstSolventAverageMass());

    double Eab = 0.0;
    if (chemBasis.snyderApproximation()) 
    {
        Eab = Nb * chemBasis.secondSolventBindEnergy();
    }
    else
    {
        Eab = 0.0 + 1.0 / columnRelativeStrength * log( 1.0 - Nb + Nb * Q );
    }

    std::vector<double> monomerEnergyProfile;
    for (std::vector<BioLCCC::ChemicalGroup>::const_iterator residue =
                ++(parsedSequence.begin());
            residue != --(parsedSequence.end());
            ++residue)
    {
        double residueEnergy = residue->bindEnergy();

        // Adding the energy of the N-terminal group to the first residue.
        if (residue == ++(parsedSequence.begin()))
        {
            residueEnergy += parsedSequence.begin()->bindEnergy();
        }

        // Adding the energy of the C-terminal group to the last residue.
        else if (residue == --(--(parsedSequence.end())))
        {
            residueEnergy += (--(parsedSequence.end()))->bindEnergy();
        }

        monomerEnergyProfile.push_back(
            columnRelativeStrength*(residueEnergy - Eab) * 293.0 / temperature);
    }
    return monomerEnergyProfile;
}

/*!
    This function calculates the effective energy profile of a polymer chain.
    Each element in this profile contains the adsorption energy of a single
    segment.

    The energy of a single segment equals to sum of the energies of the whole
    monomers containing in the segment PLUS the proportional parts of the
    energies of monomers which cross the borders of the segment. 
*/
std::vector<double> calculateSegmentEnergyProfile(
    const std::vector<double> &monomerEnergyProfile,
    const double monomerLength,
    const double kuhnLength)
{
    std::vector<double> segmentEnergyProfile; 
    std::vector<double>::const_iterator monomerEnergyIt =
        monomerEnergyProfile.begin();
    double kuhnLeftBorder  = 0;
    double monomerLeftBorder  = 0;
    double sumEnergy = 0.0;
    bool kuhnSegmentOpen = true;

    while (monomerEnergyIt != monomerEnergyProfile.end()) 
    {
        if ((kuhnLeftBorder + kuhnLength) >= 
                (monomerLeftBorder + monomerLength))
        {
            sumEnergy += *(monomerEnergyIt) * 
                (monomerLeftBorder + monomerLength - 
                    std::max(monomerLeftBorder, kuhnLeftBorder)) / 
                monomerLength;
            kuhnSegmentOpen = true;
            monomerLeftBorder += monomerLength;
            ++monomerEnergyIt;
        }
        else 
        {
            sumEnergy += *(monomerEnergyIt) * 
                (kuhnLeftBorder + kuhnLength - 
                    std::max(monomerLeftBorder, kuhnLeftBorder)) / 
                monomerLength;
            segmentEnergyProfile.push_back(sumEnergy);
            sumEnergy = 0.0;
            kuhnSegmentOpen = false;
            kuhnLeftBorder += kuhnLength;
        }
    }

    if (kuhnSegmentOpen)
    {
        segmentEnergyProfile.push_back(sumEnergy);
    }

    return segmentEnergyProfile;
}

/*!
    This function converts the energy profile of a peptide/protein to 
    a profile of distribution probabilities. Probability = exp(E_effective).
*/
std::vector<double> calculateBoltzmannFactorProfile(
    const std::vector<double> &effectiveEnergyProfile)
{
    std::vector<double> boltzmannFactorProfile; 

    for (std::vector<double>::const_iterator energy =
                effectiveEnergyProfile.begin();
            energy != effectiveEnergyProfile.end();
            ++energy)
    {
        boltzmannFactorProfile.push_back(exp(*energy));
    }

    return boltzmannFactorProfile;
}

double calculateKdChain(
    const std::vector<ChemicalGroup> &parsedSequence,
    const double secondSolventConcentration,
    const ChemicalBasis &chemBasis,
    const double columnPoreSize,
    const double columnRelativeStrength,
    const double temperature) throw (BioLCCCException)
{
    // At first, we need to convert the energy profile to a profile of 
    // distribution probabilities. Probability = exp(E_effective),
    // where E_effective = E_of_residue - Eab,
    // and Eab is an energy of binding for a binary solvent
    // and Eab = Ea + ln ( 1 + Nb + Nb*exp (Ea - Eb) )
    // also corrections of energies due to temperature (energies in exponents
    // are scaled to the new temperature) and the relative strength 
    // are introduced.
    std::vector<std::vector<double> > boltzmannFactorProfiles;
    for (std::vector<double>::const_iterator adsorptionLayerFactor = 
            chemBasis.adsorptionLayerFactors().begin();
         adsorptionLayerFactor != chemBasis.adsorptionLayerFactors().end();
         ++adsorptionLayerFactor) 
    {
        boltzmannFactorProfiles.push_back(
            calculateBoltzmannFactorProfile(
                calculateSegmentEnergyProfile(
                    calculateMonomerEnergyProfile(
                        parsedSequence,
                        chemBasis,
                        secondSolventConcentration,
                        columnRelativeStrength * (*adsorptionLayerFactor),
                        temperature),
                    chemBasis.monomerLength(),
                    chemBasis.kuhnLength())));
    }

    // The size of the lattice must be greater than 
    // (number of adsorbing layers) * 2.
    // double round (double x) {return floor(x+0.5);}
    const unsigned int latticeSize = 
        floor(columnPoreSize / chemBasis.kuhnLength() + 0.5);

    if (latticeSize < chemBasis.adsorptionLayerFactors().size() * 2)
    {
        throw BioLCCCException(
            "The pore size is too small for the given number of adsorbing "
            "layers.");
    }

    // The density vector correspond to a probability of n-th residue to be in
    // a certain layer between pore walls.
    double *density;

    // The transition matrix used to calculate a density vector of n-th
    // particle from a density vector of (n-1)-th particle.
    double *transitionMatrix;

    // The density buffer vector is used during matrix calculations.
    double *densityBuffer;

    // Memory managment.
    try
    {
        density = new double[latticeSize];
        densityBuffer = new double[latticeSize];
        transitionMatrix = new double[latticeSize*latticeSize];
    }
    catch (...)
    {
        throw BioLCCCException("Cannot allocate memory for calculations");
    }

    // Constructing a density vector for the first amino acid residue.
    // A density is distributed uniformly over all non-adsorbing layers of 
    // the lattice.
    // The density in adsorbing layers is multiplied by a Boltzmann factor of
    // the first segment.

    for (unsigned int i = 0; i < latticeSize; i++)
    {
        density[i] = 1.0;
    }

    for (unsigned int i = 0; i < boltzmannFactorProfiles.size(); ++i) 
    {
        density[i] = boltzmannFactorProfiles[i][0];
        density[latticeSize - i - 1] = boltzmannFactorProfiles[i][0];
    }

    // Debugging facilities.
    //for (unsigned int i = 0; i < latticeSize; i++)
    //{
    //    std::cout << density[i] << " ";
    //}
    //std::cout << std::endl;
    //std::cout << std::endl;

    // Than we construct a basis for the transition matrix. The basis is
    // a diagonal matrix with 4.0/6.0 on the main diagonal and 1.0/6.0 on
    // the side diagonals.

    // Filling the matrix.
    for (unsigned int i = 0; i < latticeSize; i++)
    {
        for (unsigned int j = 0; j < latticeSize; j++)
        {
            switch (abs( j - i ))
            {
            case 0:
            {
                transitionMatrix[j + latticeSize * i] = 4.0/6.0;
                break;
            }
            case 1:
            {
                transitionMatrix[j + latticeSize * i] = 1.0/6.0;
                break;
            }
            default:
                transitionMatrix[j + latticeSize * i] = 0.0;
            }
        }
    }

    // On each step we obtain the density vector for the n-th segment
    // multiplying the transition matrix and the density vector of the 
    // (n-1)th residue.
    // The multiplication starts from the second segment.
    for (unsigned int segmentIndex = 1; 
         segmentIndex < boltzmannFactorProfiles[0].size();
         ++segmentIndex) 
    {
        // Filling the matrix elements that describe the adsorption.
        for (unsigned int layerIndex = 0;
             layerIndex < boltzmannFactorProfiles.size();
             ++layerIndex) 
        {
            int indexShift = layerIndex * ( latticeSize + 1 );
            double boltzmannFactor = 
                boltzmannFactorProfiles[layerIndex][segmentIndex];

            transitionMatrix[indexShift + 0] = 4.0/6.0 * boltzmannFactor;
            transitionMatrix[indexShift + 1] = 1.0/6.0 * boltzmannFactor;
            transitionMatrix[latticeSize*latticeSize - 1 - indexShift - 0] =
                    4.0 / 6.0 * boltzmannFactor;
            transitionMatrix[latticeSize*latticeSize - 1 - indexShift - 1] =
                    1.0 / 6.0 * boltzmannFactor;

            // A segment may enter the second and further adsorption layers from
            // the inner layer (i.e. the layer lying closer to the walls).
            if (layerIndex > 0) 
            {
                transitionMatrix[indexShift - 1] = 1.0/6.0 * boltzmannFactor;
                transitionMatrix[latticeSize*latticeSize - 1 - indexShift + 1] =
                        1.0 / 6.0 * boltzmannFactor;
            }
        }

        // Zeroing the calculation buffer.
        for (unsigned int i = 0; i < latticeSize; i++)
        {
            densityBuffer[i] = 0.0;
        }

        // Multiplying the transition matrix by the density vector. The result
        // is stored in the buffer vector.
        for (unsigned int i = 0; i < latticeSize; i++)
        {
            for (unsigned int j = 0; j < latticeSize; j++)
            {
                densityBuffer[i] = densityBuffer[i] + density[j] *
                                   transitionMatrix[j + i * latticeSize];
            }
        }

        // Transferring the results from the density vector.
        for (unsigned int i = 0; i < latticeSize; i++)
        {
            density[i] = densityBuffer[i];
        }
    }

    // Kd is a sum of elements of the density vector, normalized to the size 
    // of the lattice.
    double Kd=0;
    for (unsigned int i=0; i < latticeSize; i++)
    {
        Kd += density[i];
    }
    Kd = Kd / (double)(latticeSize);

    // Cleaning memory.
    try
    {
        delete[] density;
        delete[] densityBuffer;
        delete[] transitionMatrix;
    }
    catch (...)
    {
        throw BioLCCCException("Cannot allocate memory for calculations");
    }

    return Kd;
}

double rodAdsorptionEnergy(const std::vector<double> & rodEnergyProfile,
                           unsigned int n,
                           bool reversed = false) throw(BioLCCCException)
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

double partitionFunctionRodSubmergedIntoLayer(
    double segmentLength,
    double slitWidth,
    double layerWidth,
    const std::vector<double> & rodEnergyProfile,
    bool reversed = false) throw(BioLCCCException)
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

double calculateKd(const std::vector<ChemicalGroup> &parsedSequence,
                   const double secondSolventConcentration,
                   const ChemicalBasis &chemBasis,
                   const double columnPoreSize,
                   const double columnRelativeStrength,
                   const double temperature) throw(BioLCCCException)
{
    // assymetricCalculations shows whether the Kd for the reversed molecule 
    // will differ. It happens when a molecule cannot be divided into an integer
    // number of Kuhn segments.
    bool assymetricCalculations = 
        (fmod(chemBasis.monomerLength() * parsedSequence.size(), 
              chemBasis.kuhnLength()) != 0);
    // Choosing the appropriate polymerModel.
    if (chemBasis.polymerModel()==CHAIN)
    {
        double Kd = calculateKdChain(parsedSequence,
                                    secondSolventConcentration, 
                                    chemBasis, columnPoreSize,
                                    columnRelativeStrength, temperature);
        if (assymetricCalculations) 
        {
            std::vector<ChemicalGroup> revParsedSequence = 
                parsedSequence;
            std::reverse(revParsedSequence.begin(),
                         revParsedSequence.end());
            Kd = (Kd + calculateKdChain(revParsedSequence,
                                       secondSolventConcentration,
                                       chemBasis,
                                       columnPoreSize,
                                       columnRelativeStrength, 
                                       temperature)) / 2.0 ;
        }
        return Kd;
    }
    else if (chemBasis.polymerModel() == ROD)
    {
        double Kd = calculateKdRod(parsedSequence,
                                   secondSolventConcentration, 
                                   chemBasis, columnPoreSize,
                                   columnRelativeStrength, temperature);
        if (assymetricCalculations) 
        {
            std::vector<ChemicalGroup> revParsedSequence = 
                parsedSequence;
            std::reverse(revParsedSequence.begin(),
                         revParsedSequence.end());
            Kd = (Kd + calculateKdRod(revParsedSequence,
                                      secondSolventConcentration,
                                      chemBasis,
                                      columnPoreSize,
                                      columnRelativeStrength, 
                                      temperature)) / 2.0 ;
        }
        return Kd;
    }
    else
    {
        throw BioLCCCException("Model error.");
    }
}

class KdCalculator
{
public:
    KdCalculator(const std::vector<ChemicalGroup> &parsedSequence,
                 const ChemicalBasis &chemBasis, 
                 const double columnPoreSize,
                 const double columnRelativeStrength,
                 const double temperature,
                 const int numInterpolationPoints) :
                     mParsedSequence(parsedSequence),
                     mChemicalBasis(chemBasis),
                     mColumnPoreSize(columnPoreSize),
                     mColumnRelativeStrength(columnRelativeStrength),
                     mTemperature(temperature),
                     mN(numInterpolationPoints)
    {
        if (mN > 0)
        {
            mSecondSolventConcentrations = new double[mN];
            mLogKds = new double[mN];
            // The number of extra points in the terminal segments.
            // This points significantly increase the accuracy of spline
            // interpolation.
            int NETP = 1;
            for (int i=0; i < mN; i++) 
            {
                if (i <= NETP)
                {
                    mSecondSolventConcentrations[i] =
                        i * 100.0 / (mN - 2.0 * NETP - 1.0) / (NETP + 1.0);
                }
                else if (i > (mN - NETP - 2))
                {
                    mSecondSolventConcentrations[i] = 
                        ((mN - 2.0 * NETP - 2.0)
                            + (i - mN + NETP + 2.0) / (NETP + 1.0))
                         * 100.0 / (mN - 2.0 * NETP - 1.0);
                }
                else
                {
                    mSecondSolventConcentrations[i] =
                        (i - NETP) * 100.0 / (mN - 2.0 * NETP - 1.0);
                }
                mLogKds[i] = log(calculateKd(mParsedSequence,
                    mSecondSolventConcentrations[i],
                    mChemicalBasis, 
                    mColumnPoreSize,
                    mColumnRelativeStrength,
                    mTemperature));
            }

            mSecondDers = new double[mN];
            fitSpline(mSecondSolventConcentrations, mLogKds,
                mN, mSecondDers);
        }
    }

    ~KdCalculator()
    {
        if (mN > 0)
        {
            delete[] mSecondSolventConcentrations;
            delete[] mLogKds;
            delete[] mSecondDers;
        }
    }

    double operator()(double secondSolventConcentration) 
        throw (BioLCCCException)
    {
        if (mN == 0) 
        {
            return calculateKd(mParsedSequence,
                               secondSolventConcentration,
                               mChemicalBasis, 
                               mColumnPoreSize,
                               mColumnRelativeStrength,
                               mTemperature);
        }
        else 
        {
            return exp(calculateSpline(mSecondSolventConcentrations, mLogKds,
                mSecondDers, mN, 
                secondSolventConcentration));
            //return exp(partPolInterpolate(
            //    mSecondSolventConcentrations, mLogKds,
            //    mN, 2,
            //    secondSolventConcentration));
        }
    }

private:
    const std::vector<ChemicalGroup> &mParsedSequence;
    const ChemicalBasis & mChemicalBasis;
    const double mColumnPoreSize;
    const double mColumnRelativeStrength;
    const double mTemperature;
    const int mN;
    double * mSecondSolventConcentrations;
    double * mLogKds;
    double * mSecondDers;
};

double calculateRT(const std::vector<ChemicalGroup> &parsedSequence,
                   const ChemicalBasis &chemBasis,
                   const ChromoConditions &conditions,
                   const int numInterpolationPoints,
                   const bool continueGradient,
                   const bool backwardCompatibility) throw(BioLCCCException)
{
    // Calculating column volumes.
    if (numInterpolationPoints < 0)
    {
        throw BioLCCCException(
            "The number of interpolation points must be non-negative.");
    }
    double volumeLiquidPhase = conditions.columnDiameter() *
                               conditions.columnDiameter() * 3.1415 * 
                               conditions.columnLength() / 4.0 /
                               1000.0 * (conditions.columnPorosity()-
                                   conditions.columnVpToVtot());
    double volumePore = conditions.columnDiameter() *
                        conditions.columnDiameter() * 3.1415 * 
                        conditions.columnLength() / 4.0 /
                        1000.0 * conditions.columnVpToVtot();

    // Recalculating dV. By default dV is calculated as the flow rate divided
    // by 20.
    double dV;
    if (conditions.dV()!= 0.0)
    {
        dV = conditions.dV();
    }
    else
    {
        dV = conditions.flowRate() / 20.0;
    }

    // A gradient should contain at least two points.
    if (conditions.gradient().size() < 2)
    {
        throw BioLCCCException(
            "The gradient must contain at least two points.");
    }

    // Converting the x-coordinate of the gradient to the scale of iterations
    // and the y-coordinate to the scale of second solvent concentration.
    std::vector<std::pair<int, double> > convertedGradient;

    for (Gradient::size_type i = 0;
        i != conditions.gradient().size();
        i++)
    {
        convertedGradient.push_back(
            std::pair<int,double>(
                int(floor(conditions.gradient()[i].time() *
                          conditions.flowRate() / dV)),
                (100.0 - conditions.gradient()[i].concentrationB())/100.0 *
                conditions.secondSolventConcentrationA() +
                conditions.gradient()[i].concentrationB() / 100.0 *
                conditions.secondSolventConcentrationB()));
    }

    std::vector<std::pair<int, double> >::const_iterator currentGradientPoint=
        convertedGradient.begin();
    std::vector<std::pair<int, double> >::const_iterator previousGradientPoint=
        convertedGradient.begin();

    KdCalculator kdCalculator(parsedSequence, chemBasis,
                              conditions.columnPoreSize(),
                              conditions.columnRelativeStrength(),
                              conditions.temperature(),
                              numInterpolationPoints);

    double secondSolventConcentration = 0.0;
    // The part of a column passed by molecules. When it exceeds 1.0,
    // molecule elute from the column.
    double S = 0.0;
    double dS = 0.0;
    // The current iteration number.
    int j = 0;
    while (S < 1.0)
    {
        j++;
        if (j > currentGradientPoint->first)
        {
            if (currentGradientPoint != --convertedGradient.end())
            {
                previousGradientPoint = currentGradientPoint;
                ++currentGradientPoint;

                // We could calculate the isocratic part of gradient more
                // efficiently due to the constancy of Kd.
                if (currentGradientPoint->second ==
                        previousGradientPoint->second)
                {
                    // One case is that a peptide elutes during this section or
                    // the section is the last. Then we can calculate the
                    // retention time directly.
                    bool peptideElutes = 
                        ((1.0 - S) / dV
                            * kdCalculator(currentGradientPoint->second)
                            * volumePore < 
                        (currentGradientPoint->first - j + 1));

                    if (peptideElutes ||
                        (currentGradientPoint == --convertedGradient.end()))
                    {
                        dS = dV / kdCalculator(currentGradientPoint->second) 
                            / volumePore;
                        j += (int) ceil((1.0 - S) / dS);
                        S += (int) ceil((1.0 - S) / dS) * dS;
                        break;
                    }

                    // Another case is that this section is not long enough for
                    // a peptide to elute. In this case we calculate the
                    // increase of S during this section.
                    else
                    {
                        S += dV / kdCalculator(currentGradientPoint->second)
                             / volumePore * (currentGradientPoint->first -
                                             previousGradientPoint->first);
                        j = currentGradientPoint->first;
                    }
                }
            }
            // If j exceeds the last point of a gradient, the value of the
            // second solvent concentration is calculated by a prolongation of
            // the last gradient section.
            else if (!continueGradient)
            {
                break;
            }
        }

        secondSolventConcentration = currentGradientPoint->second -
            (currentGradientPoint->second - previousGradientPoint->second) /
            (currentGradientPoint->first - previousGradientPoint->first) *
            (currentGradientPoint->first - (double) (j-0.5));
        dS = dV / kdCalculator(secondSolventConcentration) / volumePore;
        S += dS;
    }

    double RT = j * dV / conditions.flowRate() + conditions.delayTime() +
                volumeLiquidPhase / conditions.flowRate();
    if (!backwardCompatibility)
    {
        RT -= (S - 1.0) / dS * dV / conditions.flowRate();
    }
    return RT;
}
}

std::vector<ChemicalGroup> parseSequence(
    const std::string &source,
    const ChemicalBasis &chemBasis) throw(BioLCCCException)
{
    std::vector<ChemicalGroup> parsedSequence;
    ChemicalGroup NTerminus;
    ChemicalGroup CTerminus;

    // At first we need to strip the sequence from adjacent amino acids.
    // If a source sequence contains them, it should contain two dots, so we
    // need the part of sequence between them.
    std::size_t firstDotPosition = 0;
    std::size_t secondDotPosition = 0;

    // We'll use this variable to contain peptide sequence without adjacent
    // amino acids.
    std::string strippedSource = source;

    firstDotPosition = source.find(".");

    if (firstDotPosition != std::string::npos)
    {
        secondDotPosition = source.find(".", firstDotPosition+1);
        if (secondDotPosition != std::string::npos)
        {

            // If a source sequence contains more that two dots, it's broken.
            if (source.find(".", secondDotPosition+1) != std::string::npos)
            {
                throw ParsingException(
                    "The sequence " + source +" contains more than two dots.");
            }
            else
            {
                strippedSource = source.substr(firstDotPosition+1,
                    secondDotPosition - firstDotPosition - 1);
            }
        }
        // If a source sequence contains only one dot, it's broken.
        else
        {
            throw ParsingException(
                "The sequence " + source + " contains only one dot.");
        }
    }

    // Than goes parsing.
    std::size_t NTerminusPosition = 0;

    // First we need to check the stripped source sequence for
    // the N-Terminal group.
    NTerminus = chemBasis.defaultNTerminus();
    for (std::map<std::string,ChemicalGroup>::const_iterator
            NTerminusIterator = chemBasis.chemicalGroups().begin();
            NTerminusIterator != chemBasis.chemicalGroups().end();
            NTerminusIterator++)
    {
        if (NTerminusIterator->second.isNTerminal())
        {
            if (strippedSource.find(NTerminusIterator->second.label()) ==
                    (size_t)0)
            {
                NTerminus = NTerminusIterator->second;
                NTerminusPosition = NTerminusIterator->second.label().size();
            }
        }
    }

    // Then we need to found the location of the C-Terminus.
    CTerminus = chemBasis.defaultCTerminus();
    std::size_t CTerminusPosition;
    CTerminusPosition = strippedSource.find("-", NTerminusPosition);
    if (CTerminusPosition != std::string::npos)
    {
        // The sequence should not contain hyphens after C-terminal group.
        if (strippedSource.find("-", CTerminusPosition+1) != std::string::npos)
        {
            throw ParsingException(
                "The sequence " + source +
                " contains hyphen after C-terminal group.");
        }

        // Searching for known C-terminal groups.
        bool CTerminusFound = false;
        for (std::map<std::string,ChemicalGroup>::const_iterator
                CTerminusIterator = chemBasis.chemicalGroups().begin();
                CTerminusIterator != chemBasis.chemicalGroups().end();
                CTerminusIterator++)
        {
            if (CTerminusIterator->second.isCTerminal())
            {
                if (strippedSource.find(CTerminusIterator->second.label(),
                    CTerminusPosition) != std::string::npos)
                {
                    CTerminus = CTerminusIterator->second;
                    CTerminusFound = true;
                }
            }
        }
        if (!CTerminusFound)
        {
            throw ParsingException(
                "The sequence " + source +
                " contains unknown C-terminal group\"" + 
                source.substr(CTerminusPosition) + "\".");
        }
    }
    else
    {
        CTerminusPosition = strippedSource.size();
    }

    // At this step we obtain the sequence of a peptide without adjacent
    // amino acids and terminal groups.
    strippedSource = strippedSource.substr(
        NTerminusPosition, CTerminusPosition-NTerminusPosition);

    // We need to check whether it contains any non-letter characters.
    for (std::size_t i=0; i<strippedSource.size(); i++)
    {
        if (!(((int(strippedSource[i]) >= int('a')) &&
                (int(strippedSource[i]) <= int('z'))) ||
                ((int(strippedSource[i]) >= int('A')) &&
                 (int(strippedSource[i]) <= int('Z')))))
        {
            throw ParsingException(
                "The sequence " + source +
                " contains a non-letter character.");
        }
    }

    // Then we divide the whole sequence into aminoacids.
    bool aminoAcidFound;
    size_t curPos = 0;
    while (curPos < strippedSource.size())
    {
        aminoAcidFound = false;
        for (std::map<std::string,ChemicalGroup>::const_iterator
                aminoAcidIterator = chemBasis.chemicalGroups().begin();
                aminoAcidIterator != chemBasis.chemicalGroups().end();
                aminoAcidIterator++)
        {
            if (strippedSource.compare(curPos,
                aminoAcidIterator->second.label().size(),
                aminoAcidIterator->second.label()) == 0)
            {
                curPos += aminoAcidIterator->second.label().size();
                parsedSequence.push_back(aminoAcidIterator->second);
                aminoAcidFound = true;
                break;
            }
        }

        if (!aminoAcidFound)
        {
            throw ParsingException(
                "The sequence " + source + " contains unknown amino acid \"" + 
                source.substr(curPos, 1) + "\".");
        }
    }
    parsedSequence.insert(parsedSequence.begin(), NTerminus);
    parsedSequence.push_back(CTerminus);
    return parsedSequence;

}

double calculateRT(const std::string &sequence,
                   const ChemicalBasis &chemBasis,
                   const ChromoConditions &conditions,
                   const int numInterpolationPoints,
                   const bool continueGradient,
                   const bool backwardCompatibility) 
                   throw(BioLCCCException)
{
    std::vector<ChemicalGroup> parsedSequence = 
        parseSequence(sequence, chemBasis);
    return calculateRT(parsedSequence,
                       chemBasis,
                       conditions,
                       numInterpolationPoints,
                       continueGradient,
                       backwardCompatibility);
}

double calculateKd(const std::string &sequence,
                   const double secondSolventConcentration,
                   const ChemicalBasis & chemBasis,
                   const double columnPoreSize,
                   const double columnRelativeStrength,
                   const double temperature)
                   throw(BioLCCCException)
{
    return calculateKd(parseSequence(sequence, chemBasis),
                       secondSolventConcentration,
                       chemBasis,
                       columnPoreSize,
                       columnRelativeStrength,
                       temperature);
}

double calculateAverageMass(const std::string &sequence,
                            const ChemicalBasis &chemBasis)
                            throw(BioLCCCException)
{
    std::vector<ChemicalGroup> parsedSequence =
        parseSequence(sequence, chemBasis);
    double peptideAverageMass = 0;
    for (std::vector<ChemicalGroup>::const_iterator i =
                parsedSequence.begin();
            i < parsedSequence.end();
            i++)
    {
        peptideAverageMass += i -> averageMass();
    }

    return peptideAverageMass;
}

double calculateMonoisotopicMass(const std::string &sequence,
                                 const ChemicalBasis &chemBasis)
                                 throw(BioLCCCException)
{
    std::vector<ChemicalGroup> parsedSequence =
        parseSequence(sequence, chemBasis);
    double peptideMonoisotopicMass = 0;
    for (std::vector<ChemicalGroup>::const_iterator i =
                parsedSequence.begin();
            i < parsedSequence.end();
            i++)
    {
        peptideMonoisotopicMass += i -> monoisotopicMass();
    }

    return peptideMonoisotopicMass;
}

}
