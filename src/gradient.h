#ifndef GRADIENT_H
#define GRADIENT_H

#include <vector>
#include "gradientpoint.h" 

//! The Gradient class encapsulates all properties of an elution gradient.

namespace BioLCCC {

class Gradient: public std::vector<GradientPoint> {
    public:
        /*!
            Constructs an empty elution gradient.
        */
        Gradient();

        /*!
            Constructs linear elution gradient.
        */
        Gradient(double initialConcentrationB,
                 double finalConcentrationB,
                 double time);

        /*!
            Adds a new point to the gradient.
        */
        Gradient addPoint(GradientPoint iPoint);

        /*!
            Adds a new point to the gradient in more convenient way.
        */
        Gradient addPoint(double iTime, double iConcentrationB);

    private:
        double mTime;
        double mConcentrationB;
};

}

#endif
