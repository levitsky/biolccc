#include "gradient.h" 

Gradient::Gradient() {
}

Gradient Gradient::addPoint(GradientPoint iPoint) {
    //TODO: add exceptions here!!
    // Each new point should be later than the previous one.
    if ((this->size() == 0) || (iPoint.time() >= this->back().time())) {
        this->push_back(iPoint);
    }

    return (*this);
}

Gradient Gradient::addPoint(double iTime, double iConcentrationB) {
    this->addPoint(GradientPoint(iTime, iConcentrationB));
    return (*this);
}

