#include "chromoconditions.h"

namespace BioLCCC {
ChromoConditionsException::ChromoConditionsException(std::string message):
    BioLCCCException(message) {};

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
    setColumnLength(iColumnLength);
    setColumnDiameter(iColumnDiameter);
    setColumnPoreSize(iColumnPoreSize);
    setGradient(iGradient);
    setColumnVpToVtot(iColumnVpToVtot);
    setColumnPorosity(iColumnPorosity);
    setTemperature(iTemperature);
    setCalibrationParameter(iCalibrationParameter);
    setFlowRate(iFlowRate);
    setDV(iDV);
    setDelayTime(iDelayTime);
    setSecondSolventConcentrationA(iSecondSolventConcentrationA);
    setSecondSolventConcentrationB(iSecondSolventConcentrationB);
}
 
double ChromoConditions::columnLength() const {
    return mColumnLength;
}

void ChromoConditions::setColumnLength(double newColumnLength) {
    if (newColumnLength < 0.0) {
        throw(ChromoConditionsException("The new column length is negative."));
    }
    mColumnLength = newColumnLength;
}

double ChromoConditions::columnDiameter() const {
    return mColumnDiameter;
}

void ChromoConditions::setColumnDiameter(double newColumnDiameter) {
    if (newColumnDiameter < 0.0) {
        throw(ChromoConditionsException("The new column diameter is negative."));
    }
    mColumnDiameter = newColumnDiameter;
}

double ChromoConditions::columnPoreSize() const {
    return mColumnPoreSize;
}

void ChromoConditions::setColumnPoreSize(double newColumnPoreSize) {
    if (newColumnPoreSize < 0.0) {
        throw(ChromoConditionsException(
            "The new column pore size is negative."));
    }
    mColumnPoreSize = newColumnPoreSize;
}

double ChromoConditions::columnVpToVtot() const {
    return mColumnVpToVtot;
}

void ChromoConditions::setColumnVpToVtot(double newColumnVpToVtot) {
    if (newColumnVpToVtot < 0.0) {
        throw(ChromoConditionsException("The new column VpToVtot is negative."));
    }
    mColumnVpToVtot = newColumnVpToVtot;
}

double ChromoConditions::columnPorosity() const {
    return mColumnPorosity;
}

void ChromoConditions::setColumnPorosity(double newColumnPorosity) {
    if (newColumnPorosity < 0.0) {
        throw(ChromoConditionsException("The new column porosity is negative."));
    }
    if (newColumnPorosity > 1.0) {
        throw(ChromoConditionsException(
            "The new column porosity is greater than 1.0"));
    }
    mColumnPorosity = newColumnPorosity;
}

double ChromoConditions::temperature() const {
    return mTemperature;
}

void ChromoConditions::setTemperature(double newTemperature) {
    if (newTemperature < 0.0) {
        throw(ChromoConditionsException("The new temperature is negative."));
    }
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
    if (newFlowRate < 0.0) {
        throw(ChromoConditionsException("The new flow rate is negative."));
    }
    mFlowRate = newFlowRate;
}

double ChromoConditions::dV() const {
    return mDV;
}

void ChromoConditions::setDV(double newDV) {
    if (newDV < 0.0) {
        throw(ChromoConditionsException("The new dV is negative."));
    }
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
    if (newSecondSolventConcentrationA < 0.0) {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A
            is negative."));
    }
    if (newSecondSolventConcentrationA > 100.0) {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A
            is greater than 100%."));
    }
    mSecondSolventConcentrationA = newSecondSolventConcentrationA;
}

double ChromoConditions::secondSolventConcentrationB() const {
    return mSecondSolventConcentrationB;
}

void ChromoConditions::setSecondSolventConcentrationB(
    double newSecondSolventConcentrationB
) {
    if (newSecondSolventConcentrationB < 0.0) {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A
            is negative."));
    }
    if (newSecondSolventConcentrationB > 100.0) {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A
            is greater than 100%."));
    }
    mSecondSolventConcentrationB = newSecondSolventConcentrationB;
}

Gradient ChromoConditions::gradient() const {
    return mGradient;
}

void ChromoConditions::setGradient(Gradient newGradient) {
    if (newGradient.size() < 2) {
        throw ChromoConditionsException(
            "The gradient must contain at least two points.");
    }
    mGradient = newGradient;
}
}

