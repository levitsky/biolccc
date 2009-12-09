#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <vector>
#include <iostream>
#include "math.h"

namespace BioLCCC {

    template <typename T>
    T vectorNorm (std::vector<T> vec) {
        T norm2 = 0.0;
        for(unsigned int i = 0; i < vec.size(); i++) {
            norm2 += (vec[i] * vec[i]);
        }
        return sqrt(norm2);
    };

    template <typename T>
    std::vector<T> VectorSum(std::vector<T> v1, std::vector<T> v2) {
        if(v1.size() != v2.size()) {
            std::cout << "VectorSum error: vector sizes are not equal.\n";
            return v1;
        }
        std::vector<T> sum(v1.size());
        for(unsigned int i = 0; i < v1.size(); i++) {
            sum[i] = v1[i]+v2[i];
        }
        return sum;
    }

    template <typename T>
    std::vector<T> VectorDiff(std::vector<T> v1, std::vector<T> v2) {
        if(v1.size() != v2.size()) {
            std::cout << "VectorDiff error: vector sizes are not equal.\n";
            return v1;
        }
        std::vector<T> diff(v1.size());
        for(unsigned int i = 0; i < v1.size(); i++) {
            diff[i] = v1[i]-v2[i];
        }
        return diff;
    }

    template<class optimizedFunctionType, class setterFunctionType>
    std::vector<double> calculateGradient(
            optimizedFunctionType optimizedFunction,
            std::vector<setterFunctionType> coordinateSetters,
            std::vector<double> point,
            std::vector<double> steps) {

        //dimension of vector parameters
        unsigned int dim = coordinateSetters.size(); 

        if ((dim != point.size()) || (dim != steps.size())) {
             std::cout << 
                "calculateGradient error: vector sizes are not equal.\n";
             return point;
        }
        std::vector<double> gradient(dim);

        double f1, f2;
        for (unsigned int i = 0; i < dim; i++)
            coordinateSetters[i](point[i]);

        for (unsigned int i = 0; i < dim; i++) {
            coordinateSetters[i](point[i]-steps[i]);
            f1 = optimizedFunction();
            coordinateSetters[i](point[i]+steps[i]);
            f2 = optimizedFunction();
            coordinateSetters[i](point[i]);
            gradient[i] = (f1 - f2)/(2. * steps[i]);
        }
        return gradient;
    };

    template<class optimizedFunctionType, class setterFunctionType>
    std::vector<double> findBestPointOnLine(
            optimizedFunctionType optimizedFunction,
            std::vector<setterFunctionType> coordinateSetters,
            std::vector<double> initialPoint,
            std::vector<double> direction,
            double firstStep,
            double epsilon) {

        //dimension of vector parameters
        unsigned int dim = coordinateSetters.size(); 

        if ((dim != initialPoint.size()) || (dim != direction.size())) {
             std::cout << 
                "findBestPointOnLine error: vector sizes are not equal.\n";
             return initialPoint;
        }

        double step = firstStep;
        double valueAtCurrentPoint, valueAtNewPoint;
        std::vector<double> currentPoint = initialPoint;

        while (step >= epsilon) {
            // Evaluating the function at current point.
            for(unsigned int i = 0; i < dim; i++) {
                coordinateSetters[i](currentPoint[i]);
            }
            valueAtCurrentPoint = optimizedFunction();  

            // Evaluating the function one step forward.
            for(unsigned int i = 0; i < dim; i++) {
                coordinateSetters[i](currentPoint[i] + step*direction[i]); 
            }
            valueAtNewPoint = optimizedFunction();

            // If the 'further' value is better...
            if(valueAtNewPoint < valueAtCurrentPoint) {
                //... then make a step forward ...
                for(unsigned int i = 0; i < dim; i++) {
                    currentPoint[i] += step*direction[i];   
                }
                step *= 2.;                                 
                //... and double the step length.
            }
            else {  
                // If not, look at the value 1 step back.
                for(unsigned int i = 0; i < dim; i++) {
                    coordinateSetters[i](currentPoint[i] - step*direction[i]);
                }
                valueAtNewPoint = optimizedFunction();
                // If the "back" value is better...
                if(valueAtNewPoint < valueAtCurrentPoint) {     
                    //...then make a step back ...
                    for(unsigned int i = 0; i < dim; i++)
                        currentPoint[i] -= step*direction[i];   
                    //...and double the step length.
                    step *= 2.; 
                }
                // Finally, if neither direction is OK,
                // decrease the step length 2-fold
                else {                             
                    step /= 2.;                                 
                }
            }
        }
        return currentPoint;
    };

    /*!
        Finds the minimum point of the optimized function with 
        the brute force method.
    */
    template<class optimizedFunctionType, class setterFunctionType>
    std::vector<double> findMinimumBruteForce (
        optimizedFunctionType optimizedFunction,
        std::vector<setterFunctionType> coordinateSetters,
        std::vector<double> lowerBounds,
        std::vector<double> upperBounds,
        std::vector<double> steps) {

        //dimension of vector parameters
        unsigned int dim = coordinateSetters.size(); 

        if ((dim != lowerBounds.size()) || 
            (dim != upperBounds.size()) ||
            (dim != steps.size())) {
            std::cout << 
                "findMinimumBruteForce error: vector sizes are not equal.\n";
            return lowerBounds;
        }

        for(unsigned int i = 0; i < dim; i++) {
            coordinateSetters[i](lowerBounds[i]);
        }
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

        //double minPoint = lowerBounds.back();
        coordinateSetters.back()(lowerBounds.back());
        std::vector<double> subMin;
        std::vector<double> minimumPoint = lowerBounds;
        double curValue = minValue;
        for(double curPoint = lowerBounds.back();
            curPoint <= upperBounds.back();
            curPoint += steps.back()) {

                // find minimum for coordinates 0 .. dim-1 at the current value of the last coordinate
                coordinateSetters.back()(curPoint);
                subMin = findMinimumBruteForce(optimizedFunction,
                                               newCoordinateSetters,
                                               newLowerBounds,
                                               newUpperBounds,
                                               newSteps);
                //std::cout << subMin.size();
                for(unsigned int i = 0; i < newCoordinateSetters.size(); i++) {
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
    };


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

        //dimension of vector parameters
        unsigned int dim = coordinateSetters.size(), count = 0; 

        if ((dim != initialPoint.size()) || (dim != steps.size())) {
             std::cout << "findMinimumGradientDescent error: vector sizes are not equal.\n";
             return initialPoint;
        }

        std::vector<double> currentPoint = initialPoint, shift(dim), gradient, currentBestPoint;
        double firstStep;
        while(true) {
            gradient = calculateGradient(optimizedFunction, coordinateSetters, currentPoint, steps);
            firstStep = vectorNorm(steps) / vectorNorm(gradient);
            currentBestPoint = findBestPointOnLine(optimizedFunction,
                                                   coordinateSetters,
                                                   currentPoint,
                                                   gradient,
                                                   firstStep,
                                                   epsilon);
            for(unsigned int i = 0; i < dim; i++)
                shift[i] = currentPoint[i] - currentBestPoint[i];

            if(vectorNorm(shift) < epsilon)
                break;
            currentPoint = currentBestPoint;
            count++;
        }
        std::cout << "findMinimumGradientDescent: count = " << count << "\n";
        return currentPoint;
    };

}

#endif
