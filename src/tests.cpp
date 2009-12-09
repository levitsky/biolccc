#include <iostream>
#include <vector>
#include "boost/function.hpp"
#include "boost/bind.hpp"
#include "BioLCCC.h"
#include "tests.h"


BioLCCC::BruteForceTester::BruteForceTester () {
    x = 0.0;
    y = 0.0;
    z = 0.0;
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
/*    BioLCCC::BruteForceTester tester;
    boost::function<double()> calc = boost::bind(&BioLCCC::BruteForceTester::calculate, &tester);
    std::cout << "Test: " << calc() << "\n";

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
    min = BioLCCC::findMinimumGradientDescent(calc, setters, initialPoint, steps, 1.0e-10);
    std::cout << "GradientDescent: Minimum at (" << min[0] << ", " << min[1] << ", " << min[2] << ").\n";
*/
    std::vector<std::string> peptides,groups;
    std::vector<double> rTimes;

    peptides.push_back("Ac-GVGKGGVGVK-NH2"); rTimes.push_back(15.48);
    peptides.push_back("Ac-VVKGGVGKVGV-NH2"); rTimes.push_back(26.13);
    peptides.push_back("Ac-KGVGKVGGVK-NH2"); rTimes.push_back(13.98);
    peptides.push_back("Ac-VVGVKGGVGK-NH2"); rTimes.push_back(21.95);
    peptides.push_back("Ac-KGVGGKVGVV-NH2"); rTimes.push_back(21.05);
    peptides.push_back("Ac-GVGGVK-NH2"); rTimes.push_back(8.58);
    peptides.push_back("Ac-KGVVKGVGVKGGVKG-NH2"); rTimes.push_back(23.88);
    groups.push_back("Ac-"); groups.push_back("G"); groups.push_back("V");
    groups.push_back("K"); groups.push_back("-NH2"); groups.push_back("ACN");

    BioLCCC::ChromoConditions chr;
    BioLCCC::ChemicalBasis chb;
    BioLCCC::ChemicalBasis calibrated = BioLCCC::calibrateBioLCCC(peptides, rTimes, chr, chb, groups);

    return 0;
}
