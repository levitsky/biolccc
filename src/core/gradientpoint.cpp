#include "gradientpoint.h"

namespace BioLCCC {

GradientPoint::GradientPoint(double time,
    double concentrationB
) {
    mTime = 0.0;
    mConcentrationB = 0.0;
    if ((time >=0) && (concentrationB >=0) && (concentrationB <=100)) {
       mTime = time;
       mConcentrationB = concentrationB;
    }
}

double GradientPoint::time() const{
    return mTime;
}

double GradientPoint::concentrationB() const{
    return mConcentrationB;
}

}

