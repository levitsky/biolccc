#include "chromoconditions.h"

// ChromoConditions::ChromoConditions() {
//      // Setting a protocol "Pilot".
//     mColumnLength = 150;
//     mColumnDiameter = 0.075;
//     mColumnPoreSize = 100;
//     mColumnVpToVtot = 0.5;
//     mColumnPorosity = 0.9;
//     mTemperature = 293;
//     mCalibrationParameter = 1;
//     mFlowRate = 0.0003;
//     mDV = 0;
//     mDelayTime = 0; 
//     mSecondSolventConcentrationA = 2;
//     mSecondSolventConcentrationB = 80;
//     mSecondSolvent = std::string("ACN");
//     mGradient.push_back(gradientPoint(0,0));
//     mGradient.push_back(gradientPoint(60,50));
// }

ChromoConditions::ChromoConditions(double iColumnLength,
                                double iColumnDiameter,
                                double iColumnPoreSize,
                                gradientProfile iGradient,
                                double iSecondSolventConcentrationA,
                                double iSecondSolventConcentrationB,
                                double iDelayTime,
                                std::string iSecondSolvent,
                                double iFlowRate,
                                double iDV,
                                double iCalibrationParameter,
                                double iColumnVpToVtot,
                                double iColumnPorosity,
                                double iTemperature
) {
    mColumnLength = iColumnLength;
    mColumnDiameter = iColumnDiameter;
    mColumnPoreSize = iColumnPoreSize;
    mColumnVpToVtot = iColumnVpToVtot;
    mColumnPorosity = iColumnPorosity;
    mTemperature = iTemperature;
    mCalibrationParameter = iCalibrationParameter;
    mFlowRate = iFlowRate;
    mDV = iDV;
    mDelayTime = iDelayTime; 
    mSecondSolvent = iSecondSolvent;
    mSecondSolventConcentrationA = iSecondSolventConcentrationA;
    mSecondSolventConcentrationB = iSecondSolventConcentrationB;
    
    //setting default gradient.
    if (iGradient.empty()) {
        mGradient = gradientProfile();
        mGradient.push_back(gradientPoint(0,0));
        mGradient.push_back(gradientPoint(60,50));
    }
    else {
        mGradient = iGradient;
    }
}
 
double ChromoConditions::columnLength() const {
    return mColumnLength;
}

void ChromoConditions::setColumnLength(double newColumnLength) const {
    mColumnLength = newColumnLength;
}

double ChromoConditions::columnDiameter() const {
    return mColumnDiameter;
}

void ChromoConditions::setColumnDiameter(double newColumnDiameter) const {
    mColumnDiameter = newColumnDiameter;
}

double ChromoConditions::columnPoreSize() const {
    return mColumnPoreSize;
}

void ChromoConditions::setColumnPoreSize(double newColumnPoreSize) const {
    mColumnPoreSize = newColumnPoreSize;
}

double ChromoConditions::columnVpToVtot() const {
    return mColumnVpToVtot;
}

void ChromoConditions::setColumnVpToVtot(double newColumnVpToVtot) const {
    mColumnVpToVtot = newColumnVpToVtot;
}

double ChromoConditions::columnPorosity() const {
    return mColumnPorosity;
}

void ChromoConditions::setColumnPorosity(double newColumnPorosity) const {
    mColumnPorosity = newColumnPorosity;
}

double ChromoConditions::temperature() const {
    return mTemperature;
}

void ChromoConditions::setTemperature(double newTemperature) const {
    mTemperature = newTemperature;
}

double ChromoConditions::calibrationParameter() const {
    return mCalibrationParameter;
}

void ChromoConditions::setCalibrationParameter(
                       double newCalibrationParameter) const {
    mCalibrationParameter = newCalibrationParameter;
}

double ChromoConditions::flowRate() const {
    return mFlowRate;
}

void ChromoConditions::setFlowRate(double newFlowRate) const {
    mFlowRate = newFlowRate;
}

double ChromoConditions::dV() const {
    return mDV;
}

void ChromoConditions::setDV(double newDV) const {
    mDV = newDV;
}

double ChromoConditions::delayTime() const {
    return mDelayTime;
}

void ChromoConditions::setDelayTime(double newDelayTime) const {
    mDelayTime = newDelayTime;
}

std::string ChromoConditions::secondSolvent() const {
    return mSecondSolvent;
}

double ChromoConditions::secondSolventConcentrationA() const {
    return mSecondSolventConcentrationA;
}

double ChromoConditions::secondSolventConcentrationB() const {
    return mSecondSolventConcentrationB;
}

std::vector<std::pair<double,double> >::const_iterator ChromoConditions::beginGradient() const {
    return mGradient.begin();
}

std::vector<std::pair<double,double> >::const_iterator ChromoConditions::endGradient() const {
    return mGradient.end();
}
