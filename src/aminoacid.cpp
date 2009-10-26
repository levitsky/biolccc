#include "aminoacid.h"

Aminoacid :: Aminoacid(std::string Name,
    std::string Label,
    double bindEnergy,
    double averageMass,
    double monoisotopicMass
) {
    mName = Name;
    mLabel = Label;
    mAverageMass = averageMass;
    mMonoisotopicMass = monoisotopicMass;
    mBindEnergy = bindEnergy;
}

std::string Aminoacid::name() const{
    return mName;
}

std::string Aminoacid::label() const{
    return mLabel;
}

double Aminoacid::bindEnergy() const{
    return mBindEnergy;
}

double Aminoacid::averageMass() const{
    return mAverageMass;
}

double Aminoacid::monoisotopicMass() const{
    return mMonoisotopicMass;
}
