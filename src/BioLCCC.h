#ifndef PEPTIDEMETHODS_H
#define PEPTIDEMETHODS_H

#include "chemicalbasis.h"
#include "chromoconditions.h"

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
    const ChemicalBasis & chemBasis = standardChemicalBasis);

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

/*!
    Finds the minimum point of the optimized function with brute force method.
*/
template<class optimizedFunctionType, class setterFunctionType>
std::vector<double> findMinimumBruteForce (optimizedFunctionType optimizedFunction,
                                           std::vector<setterFunctionType> coordinateSetters,
                                           std::vector<double> lowerBounds,
                                           std::vector<double> upperBounds,
                                           std::vector<double> steps) {

    int dim = coordinateSetters.size(); //dimension of vector parameters

    if (! dim == lowerBounds.size() == upperBounds.size() == steps.size()) {
         std::cout << "findMinimumBruteForce error: vector sizes are not equal.\n";
         return lowerBounds;
    }

    // Initialize minValue with the value at lower bound
    for(int i = 0; i < dim; i++) coordinateSetters[i](lowerBounds[i]);
    double minValue = optimizedFunction();

//    std::cout << "optimizeCoordinates called with dim = " << dim << ".\n";

    if (dim == 1) {
        double minPoint = lowerBounds[0];
        double curValue;
        for (double curPoint = lowerBounds[0];
        curPoint <= upperBounds[0];
        curPoint += steps[0]) {
            coordinateSetters[0](curPoint);
            curValue = optimizedFunction();
            if (curValue < minValue) {
                minValue = curValue;
                minPoint = curPoint;
            }
        }
        return std::vector<double>(1, minPoint);
    }

    std::vector<setterFunctionType> newCoordinateSetters = coordinateSetters;
    std::vector<double> newLowerBounds = lowerBounds;
    std::vector<double> newUpperBounds = upperBounds;
    std::vector<double> newSteps = steps;
    newCoordinateSetters.pop_back();
    newLowerBounds.pop_back();
    newUpperBounds.pop_back();
    newSteps.pop_back();

    double minPoint = lowerBounds.back();
    coordinateSetters.back()(lowerBounds.back());
    std::vector<double> subMin;
    std::vector<double> minimumPoint;
    double curValue;
    for(double curPoint = lowerBounds.back();
        curPoint <= upperBounds.back();
        curPoint += steps.back()) {

            /* find minimum for coordinates 0 .. dim-1 and the current value of the last coordinate */
            coordinateSetters.back()(curPoint);
            subMin = findMinimumBruteForce(optimizedFunction,
                                           newCoordinateSetters,
                                           newLowerBounds,
                                           newUpperBounds,
                                           newSteps);
            for(int i = 0; i < newCoordinateSetters.size(); i++) {
                newCoordinateSetters[i](subMin[i]);
            }
            curValue = optimizedFunction();
            if (curValue < minValue) {
                minValue = curValue;
                minimumPoint = subMin;
                minimumPoint.push_back(curPoint);
            }
    }
    return minimumPoint;
}



/*
template<class optimizedFunctionType, class setterFunctionType>
std::vector<double> findMinimumBruteForce(
        optimizedFunctionType optimizedFunction,
        std::vector<setterFunctionType> coordinateSetters,
        std::vector<double> lowerBounds,
        std::vector<double> upperBounds,
        std::vector<double> steps) {

    double f, currentMinimumValue;
    std::vector<double> currentMinimumPoint, v;

    v = currentMinimumPoint = lowerBounds;
    int dim = coordinateSetters.size();

    /* check the equality of all vectors' dimensions */
 /*   if (! dim == lowerBounds.size() == upperBounds.size() == steps.size()) {
        std::cout << "findMinimumBruteForce error: vector sizes are not equal.\n";
        return lowerBounds;
    }

    /* initialize currentMinimumValue with the value at lowerBounds */
/*    for(int i = 0; i < dim; i++) coordinateSetters[i](lowerBounds[i]);
    currentMinimumValue = optimizedFunction();

//    std::cout << "dim = " << dim << "\n";
//    std::cout << "Initialized by " << currentMinimumValue << "\n";

    /* find minimum */
 /*   return optimizeCoordinates(optimizedFunction,
                        coordinateSetters,
                        lowerBounds,
                        upperBounds,
                        steps);
    /*for(int i = 0; i < dim; i++){ //'i' is the number of the current coordinate
        for(v[i] = lowerBounds[i], coordinateSetters[i](lowerBounds[i]);
        v[i] < upperBounds[i];
        coordinateSetters[i](v[i]+= steps[i])){ //screening the i-th coordinate
            f = optimizedFunction();
            if(f < currentMinimumValue){
                currentMinimumValue = f;
                currentMinimumPoint = v;
            }
        }
    }*/
//    return currentMinimumPoint;
/*}*/
}

#endif
