#include "gradientpoint.h"

namespace BioLCCC
{
GradientPointException::GradientPointException(std::string message):
        BioLCCCException(message) {};

GradientPoint::GradientPoint(double time,
                             double concentration)
                             throw (GradientPointException)
{
    setTime(time);
    setConcentration(concentration);
}

double GradientPoint::time() const
{
    return mTime;
}

double GradientPoint::concentration() const
{
    return mConcentration;
}

void GradientPoint::setTime(double newTime) 
    throw (GradientPointException)
{
    if (newTime < 0.0)
    {
        throw GradientPointException("Time is negative.");
    }

    mTime = newTime;
}

void GradientPoint::setConcentration(double newConcentration) 
    throw (GradientPointException)
{
    if (newConcentration < 0.0)
    {
        throw GradientPointException(
            "The concentration of B component is negative.");
    }

    if (newConcentration > 100.0)
    {
        throw GradientPointException(
            "The concentration of B component is greater than 100%.");
    }

    mConcentration = newConcentration;
}
}

