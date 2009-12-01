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

namespace {

double VectorNorm (std::vector<double> vector) {
double norm = 0.0;
for(int i = 0; i < vector.size(); i++)
    norm += (vector[i] * vector[i]);
return norm;
}

template<class optimizedFunctionType, class setterFunctionType>
std::vector<double> calculateGradient(
        optimizedFunctionType optimizedFunction,
        std::vector<setterFunctionType> coordinateSetters,
        std::vector<double> point,
        std::vector<double> steps) {

    int dim = coordinateSetters.size(); //dimension of vector parameters

    if (! dim == point.size() == steps.size()) {
         std::cout << "calculateGradient error: vector sizes are not equal.\n";
         return point;
    }
    std::vector<double> gradient(dim);

    double f1, f2;
    for (int i = 0; i < dim; i++)
        coordinateSetters[i](point[i]);

    for (int i = 0; i < dim; i++) {
        coordinateSetters[i](point[i]-steps[i]);
        f1 = optimizedFunction();
        coordinateSetters[i](point[i]+steps[i]);
        f2 = optimizedFunction();
        coordinateSetters[i](point[i]);
        gradient[i] = (f1 - f2)/(2. * steps[i]);
    }
    return gradient;
}

template<class optimizedFunctionType, class setterFunctionType>
std::vector<double> findBestPointOnLine(
        optimizedFunctionType optimizedFunction,
        std::vector<setterFunctionType> coordinateSetters,
        std::vector<double> initialPoint,
        std::vector<double> direction,
        double firstStep,
        double epsilon) {

    int dim = coordinateSetters.size(); //dimension of vector parameters

    if (! dim == initialPoint.size() == direction.size()) {
         std::cout << "findBestPointOnLine error: vector sizes are not equal.\n";
         return initialPoint;
    }

    double step = firstStep;
    double valueAtCurrentPoint, valueAtNewPoint;
    std::vector<double> currentPoint = initialPoint;

    while (step >= epsilon) {

        for(int i = 0; i < dim; i++) {
            coordinateSetters[i](currentPoint[i]);
        }
        valueAtCurrentPoint = optimizedFunction();  //evaluate the function at current point

        for(int i = 0; i < dim; i++) {
            coordinateSetters[i](currentPoint[i] + step*direction[i]); //evaluate the function 1 step forward
        }
        valueAtNewPoint = optimizedFunction();

        if(valueAtNewPoint < valueAtCurrentPoint) { //compare; if the 'further value is better...
            for(int i = 0; i < dim; i++)
                currentPoint[i] += step*direction[i];   //... then make a step forward
            step *= 2.;                                 //... and double the step length
        }
        else {                                          //else
            for(int i = 0; i < dim; i++) {
            coordinateSetters[i](currentPoint[i] - step*direction[i]);  //look at the value 1 step back
            }
            valueAtNewPoint = optimizedFunction();
            if(valueAtNewPoint < valueAtCurrentPoint)                   //if the value there is better...
                for(int i = 0; i < dim; i++) {
                    currentPoint[i] -= step*direction[i];               //...then make a tep back
                    step *= 2.;                                         //...and double the step length
                }
            else                                                        //or, if neither direction is OK
                step /= 2.;                                             //decrease the step length 2-fold
        }
    }
    return currentPoint;
}

}

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

    for(int i = 0; i < dim; i++) coordinateSetters[i](lowerBounds[i]);
    double minValue = optimizedFunction();

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


/*!
  Finds the local minimum of the given function with gradient descent method.
  'steps' contains increments given to coordinates for gradient evaluation.
  'epsilon' is desired (argument) accuracy.
*/
template<class optimizedFunctionType, class setterFunctionType>
std::vector<double> findMinimumGradientDescent(
        optimizedFunctionType optimizedFunction,
        std::vector<setterFunctionType> coordinateSetters,
        std::vector<double> initialPoint,
        std::vector<double> steps,
        double epsilon) {

    int dim = coordinateSetters.size(); //dimension of vector parameters

    if (! dim == initialPoint.size() == steps.size()) {
         std::cout << "findMinimumGradientDescent error: vector sizes are not equal.\n";
         return initialPoint;
    }

    std::vector<double> currentPoint = initialPoint, shift(dim), gradient, currentBestPoint;
    double firstStep;
    while(true) {
        gradient = calculateGradient(optimizedFunction, coordinateSetters, currentPoint, steps);
        firstStep = VectorNorm(steps) / VectorNorm(gradient);
        currentBestPoint = findBestPointOnLine(optimizedFunction,
                                               coordinateSetters,
                                               currentPoint,
                                               gradient,
                                               firstStep,
                                               epsilon);
        for(int i = 0; i < dim; i++)
            shift[i] = currentPoint[i] - currentBestPoint[i];
        if(VectorNorm(shift) < epsilon)
            break;
    }
    return currentPoint;
}
}

#endif
