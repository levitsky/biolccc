#ifndef GRADIENT_H
#define GRADIENT_H

#include <vector>
#include "gradientpoint.h" 

//! The Gradient class encapsulates all properties of an elution gradient.

class Gradient: public std::vector<GradientPoint> {
    public:
        /*!
            Constructs an empty elution gradient.
        */
        Gradient();

        /*!
            Adds a new point to the gradient.
        */
        Gradient addPoint(GradientPoint newPoint);
};

#endif
