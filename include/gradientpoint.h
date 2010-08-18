#ifndef GRADIENTPOINT_H
#define GRADIENTPOINT_H

#include "biolcccexception.h"

namespace BioLCCC
{

//! This exception is raised when something goes wrong with a GradientPoint.
class GradientPointException : public BioLCCCException
{
public:
    //! Constructs an instance of GradientPointException with the given message.
    GradientPointException(std::string message);
};

//! An instance of GradientPoint keeps the properties of a point of a gradient.
/*!
    Briefly, an instance of GradientPoint is a pair of values. The first one
    is time measured in minutes from the start of the gradient. The second is 
    the concentration of component B in binary solvent.  
 */

class GradientPoint
{
public:
    //! Default constructor.
    /*!
        Constructs a point of a gradient with the given time and the
        concentration of the component B of binary solution.
        \param iTime The time of the point. iTime >= 0.0
        \param iConcentrationB The concentration of the component B at the point
        in percents. 0.0 <= iConcentrationB <= 100.0
    */
    GradientPoint(double iTime = 0.0,
                  double iConcentrationB = 0.0
                 );

    //! Returns the time of the point in minutes.
    double time() const;

    //! Returns the concentration of the component B at the point in percents.
    double concentrationB() const;

private:
    double mTime;
    double mConcentrationB;
};

}

#endif

