#include "aminoacid.h"

namespace BioLCCC {

Aminoacid::Aminoacid(std::string name,
    std::string label,
    double bindEnergy,
    double averageMass,
    double monoisotopicMass
) {
    mName = name;
    mLabel = label;
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

}

