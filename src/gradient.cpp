#include "gradient.h" 

Gradient::Gradient() {
}

Gradient Gradient::addPoint(GradientPoint newPoint) {
    //TODO: add exceptions here!!
    if (newPoint.time() >= this->back().time()) {
        this->push_back(newPoint);
    }

    return (*this);
}

