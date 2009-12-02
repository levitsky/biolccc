#include <iostream>
#include <vector>
#include "boost/function.hpp"
#include "boost/bind.hpp"
#include "BioLCCC.h"
#include "tests.h"


BioLCCC::BruteForceTester::BruteForceTester () {
    x = 0.0;
    y = 0.0;
    z = 0,0;
}

double BioLCCC::BruteForceTester::calculate (void) {
    return ((x-3)*(x-3) + (y-2)*(y-2) + (z-1)*(z-1) + 2.0*((x-4)*(x-4) + (y+1)*(y+1) + (z+3)*(z+3)));
}

void BioLCCC::BruteForceTester::set_x (double new_x) {
    x = new_x;
}

void BioLCCC::BruteForceTester::set_y (double new_y) {
    y = new_y;
}

void BioLCCC::BruteForceTester::set_z (double new_z) {
    z = new_z;
}

int main () {
    BioLCCC::BruteForceTester tester;
    boost::function<double()> calc = boost::bind(&BioLCCC::BruteForceTester::calculate, &tester);
//    std::cout << "Test: " << calc() << "\n";

    std::vector< boost::function<void(double)> > setters;
    setters.push_back(boost::bind(&BioLCCC::BruteForceTester::set_x, &tester, _1));
    setters.push_back(boost::bind(&BioLCCC::BruteForceTester::set_y, &tester, _1));
    setters.push_back(boost::bind(&BioLCCC::BruteForceTester::set_z, &tester, _1));
    std::vector<double> low;
    low.push_back(-5.0); low.push_back(-5.0); low.push_back(-5.0);
    std::vector<double> high;
    high.push_back(5.0); high.push_back(5.0); high.push_back(5.0);
    std::vector<double> steps;
    steps.push_back(.05); steps.push_back(.05); steps.push_back(.05);

    std::vector<double> min;
    min = BioLCCC::findMinimumBruteForce(calc, setters, low, high, steps);
    std::cout << "BruteForce: Minimum at (" << min[0] << ", " << min[1] << ", " << min[2] << ").\n";

    std::vector<double> initialPoint;
    initialPoint.push_back(0.0); initialPoint.push_back(0.0); initialPoint.push_back(0.0);
    min = BioLCCC::findMinimumGradientDescent(calc, setters, initialPoint, steps, 0.001);
    std::cout << "GradientDescent: Minimum at (" << min[0] << ", " << min[1] << ", " << min[2] << ").\n";
    return 1;
}
