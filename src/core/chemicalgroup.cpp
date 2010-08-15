#include "chemicalgroup.h"

namespace BioLCCC {

ChemicalGroup::ChemicalGroup(std::string name,
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

std::string ChemicalGroup::name() const{
    return mName;
}

std::string ChemicalGroup::label() const{
    return mLabel;
}

double ChemicalGroup::bindEnergy() const{
    return mBindEnergy;
}

double ChemicalGroup::averageMass() const{
    return mAverageMass;
}

double ChemicalGroup::monoisotopicMass() const{
    return mMonoisotopicMass;
}

bool ChemicalGroup::isNTerminal() const{
    return (mLabel.find("-") == (size_t)(mLabel.size()-1));
}

bool ChemicalGroup::isCTerminal() const{
    return (mLabel.find("-") == (size_t)(0));
}

bool ChemicalGroup::isAminoAcid() const{
    return (!(isCTerminal() || isNTerminal()));
}

void ChemicalGroup::setBindEnergy(double newBindEnergy){
    mBindEnergy = newBindEnergy;
}
}

