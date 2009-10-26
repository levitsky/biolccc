#include "terminus.h"

Terminus :: Terminus(std::string Name,
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

std::string Terminus::name() const{
    return mName;
}

std::string Terminus::label() const{
    return mLabel;
}

double Terminus::bindEnergy() const{
    return mBindEnergy;
}

double Terminus::averageMass() const{
    return mAverageMass;
}

double Terminus::monoisotopicMass() const{
    return mMonoisotopicMass;
}
