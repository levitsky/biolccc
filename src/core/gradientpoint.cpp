#include "gradientpoint.h"

namespace BioLCCC
{
GradientPointException::GradientPointException(std::string message):
        BioLCCCException(message) {};

GradientPoint::GradientPoint(double time,
                             double concentrationB
                            )
{
    if (time < 0.0)
    {
        throw GradientPointException("Time is negative.");
    }

    if (concentrationB < 0.0)
    {
        throw GradientPointException(
            "The concentration of B component is negative.");
    }

    if (concentrationB > 100.0)
    {
        throw GradientPointException(
            "The concentration of B component is greater than 100%.");
    }

    mTime = time;
    mConcentrationB = concentrationB;
}

double GradientPoint::time() const
{
    return mTime;
}

double GradientPoint::concentrationB() const
{
    return mConcentrationB;
}

}

