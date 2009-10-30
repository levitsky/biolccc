#include "chromoconditions.h"

namespace BioLCCC {

ChromoConditions::ChromoConditions(double iColumnLength,
                                   double iColumnDiameter,
                                   double iColumnPoreSize,
                                   Gradient iGradient,
                                   double iSecondSolventConcentrationA,
                                   double iSecondSolventConcentrationB,
                                   double iDelayTime,
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
    mGradient = iGradient;
    mColumnVpToVtot = iColumnVpToVtot;
    mColumnPorosity = iColumnPorosity;
    mTemperature = iTemperature;
    mCalibrationParameter = iCalibrationParameter;
    mFlowRate = iFlowRate;
    mDV = iDV;
    mDelayTime = iDelayTime; 
    //TODO: change the behaviour of the second solvent.
    mSecondSolvent = "nothing";
    mSecondSolventConcentrationA = iSecondSolventConcentrationA;
    mSecondSolventConcentrationB = iSecondSolventConcentrationB;
}
 
double ChromoConditions::columnLength() const {
    return mColumnLength;
}

void ChromoConditions::setColumnLength(double newColumnLength) {
    mColumnLength = newColumnLength;
}

double ChromoConditions::columnDiameter() const {
    return mColumnDiameter;
}

void ChromoConditions::setColumnDiameter(double newColumnDiameter) {
    mColumnDiameter = newColumnDiameter;
}

double ChromoConditions::columnPoreSize() const {
    return mColumnPoreSize;
}

void ChromoConditions::setColumnPoreSize(double newColumnPoreSize) {
    mColumnPoreSize = newColumnPoreSize;
}

double ChromoConditions::columnVpToVtot() const {
    return mColumnVpToVtot;
}

void ChromoConditions::setColumnVpToVtot(double newColumnVpToVtot) {
    mColumnVpToVtot = newColumnVpToVtot;
}

double ChromoConditions::columnPorosity() const {
    return mColumnPorosity;
}

void ChromoConditions::setColumnPorosity(double newColumnPorosity) {
    mColumnPorosity = newColumnPorosity;
}

double ChromoConditions::temperature() const {
    return mTemperature;
}

void ChromoConditions::setTemperature(double newTemperature) {
    mTemperature = newTemperature;
}

double ChromoConditions::calibrationParameter() const {
    return mCalibrationParameter;
}

void ChromoConditions::setCalibrationParameter(double newCalibrationParameter) {
    mCalibrationParameter = newCalibrationParameter;
}

double ChromoConditions::flowRate() const {
    return mFlowRate;
}

void ChromoConditions::setFlowRate(double newFlowRate) {
    mFlowRate = newFlowRate;
}

double ChromoConditions::dV() const {
    return mDV;
}

void ChromoConditions::setDV(double newDV) {
    mDV = newDV;
}

double ChromoConditions::delayTime() const {
    return mDelayTime;
}

void ChromoConditions::setDelayTime(double newDelayTime) {
    mDelayTime = newDelayTime;
}

std::string ChromoConditions::secondSolvent() const {
    return mSecondSolvent;
}

double ChromoConditions::secondSolventConcentrationA() const {
    return mSecondSolventConcentrationA;
}

void ChromoConditions::setSecondSolventConcentrationA(
    double newSecondSolventConcentrationA
) {
    mSecondSolventConcentrationA = newSecondSolventConcentrationA;
}

double ChromoConditions::secondSolventConcentrationB() const {
    return mSecondSolventConcentrationB;
}

void ChromoConditions::setSecondSolventConcentrationB(
    double newSecondSolventConcentrationB
) {
    mSecondSolventConcentrationB = newSecondSolventConcentrationB;
}

Gradient ChromoConditions::gradient() const {
    return mGradient;
}

void ChromoConditions::setGradient(Gradient newGradient) {
    mGradient = newGradient;
}

}

