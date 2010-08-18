#ifndef GRADIENT_H
#define GRADIENT_H

#include <vector>
#include "gradientpoint.h"
#include "biolcccexception.h"

namespace BioLCCC
{

//! This exception is raised when something goes wrong with a Gradient.
class GradientException : public BioLCCCException
{
public:
    //! Constructs an instance of GradientException with a given message.
    GradientException(std::string message);
};

//! An instance of Gradient describes an elution gradient.
/*!
    An instance of Gradient keeps the shape of an elution gradient, i.e. the
    dependence of the concentration of component B on time. The shape is
    represented by a series of points (time, concentration of B component). The
    concentration in a moment of time between any two points is linearly 
    interpolated.

    Note, that gradient must start with point at 0.0 min and contain only
    successive points.

    Gradient is a public "heir" of std::vector, so see the corresponding
    reference (http://www.cplusplus.com/reference/stl/vector/) for help on 
    other members and methods of this class.
 */

class Gradient: public std::vector<GradientPoint>
{
public:
    //! Constructs an empty elution gradient.
    Gradient();

    //! Constructs a linear elution gradient.
    /*!
        Constructs a gradient with two points: (0.0, initialConcentrationB) and
        (time, finalConcentrationB).
        \param initialConcentrationB Concentration of component B at 0.0 min.
        \param finalConcentrationB Concentration of component B at the end of
        the gradient.
        \param time The duration of the gradient.
    */
    Gradient(double initialConcentrationB,
             double finalConcentrationB,
             double time);

    //! Extends the gradient with the given point.
    /*!
        \param iPoint GradientPoint to be added to the gradient.
    */
    Gradient addPoint(GradientPoint iPoint);

    //! Extends the gradient with a point at given coordinates.
    /*!
        \param iTime The time of the point. iTime >= 0.0
        \param iConcentrationB The concentration of the component B at the point
        in percents. 0.0 <= iConcentrationB <= 100.0
    */
    Gradient addPoint(double iTime, double iConcentrationB);
};

}

#endif
