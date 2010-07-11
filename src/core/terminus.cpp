#include "terminus.h"

namespace BioLCCC {
Terminus :: Terminus(std::string name,
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

void Terminus::setBindEnergy(double newBindEnergy){
    mBindEnergy = newBindEnergy;
}
}
