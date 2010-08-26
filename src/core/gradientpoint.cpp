#include "gradientpoint.h"

namespace BioLCCC
{
GradientPointException::GradientPointException(std::string message):
        BioLCCCException(message) {};

GradientPoint::GradientPoint(double time,
                             double concentrationB
                            )
{

    setTime(time);
    setConcentrationB(concentrationB);
}

double GradientPoint::time() const
{
    return mTime;
}

double GradientPoint::concentrationB() const
{
    return mConcentrationB;
}

void GradientPoint::setTime(double newTime) 
{
    if (newTime < 0.0)
    {
        throw GradientPointException("Time is negative.");
    }

    mTime = newTime;
}

void GradientPoint::setConcentrationB(double newConcentrationB) 
{
    if (newConcentrationB < 0.0)
    {
        throw GradientPointException(
            "The concentration of B component is negative.");
    }

    if (newConcentrationB > 100.0)
    {
        throw GradientPointException(
            "The concentration of B component is greater than 100%.");
    }

    mConcentrationB = newConcentrationB;
}
}

