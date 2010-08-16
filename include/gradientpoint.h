#ifndef GRADIENTPOINT_H
#define GRADIENTPOINT_H

#include "biolcccexception.h"

//! The GradientPoint class encapsulates all properties of a point of a
//gradient.

namespace BioLCCC {

class GradientPointException : public BioLCCCException {
    public:
        GradientPointException(std::string message);
};

class GradientPoint {
    public:
        /*!
            Constructs a point of a gradient with the given time and the
            concentration of the solvent B.
        */
        GradientPoint(double iTime = 0.0,
                      double iConcentrationB = 0.0
        );

        /*
            Returns the time of the point in minutes.
        */
        double time() const;

        /*
            Returns the corresponding concentration of the component B 
            in percents.
        */
        double concentrationB() const;

    private:
        double mTime;
        double mConcentrationB;
};

}

#endif

