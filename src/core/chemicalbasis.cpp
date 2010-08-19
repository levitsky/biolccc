#include "chemicalbasis.h"

namespace BioLCCC
{
ChemicalBasisException::ChemicalBasisException(std::string message):
        BioLCCCException(message) {};

ChemicalBasis::ChemicalBasis()
{
    setModel(COIL_BOLTZMANN);
    setSecondSolventBindEnergy(1.0e-10);
    setSegmentLength(1.0e-10);
    setAdsorbtionLayerWidth(1.0e-10);
    setKuhnLength(1);
}

ChemicalBasis::ChemicalBasis(PredefinedChemicalBasis predefinedChemicalBasisId)
{
    setPredefinedChemicalBasis(predefinedChemicalBasisId);
}

const std::map<std::string,ChemicalGroup> & 
    ChemicalBasis::chemicalGroups() const
{
    return mChemicalGroups;
}

const ChemicalGroup & ChemicalBasis::defaultNTerminus() const
{
    std::map<std::string, ChemicalGroup>::const_iterator NTerminusIterator =
        mChemicalGroups.find("H-");
    if (NTerminusIterator == mChemicalGroups.end())
    {
        throw ChemicalBasisException("The default H- N-terminus not found.");
    }
    return NTerminusIterator->second;
}

const ChemicalGroup & ChemicalBasis::defaultCTerminus() const
{
    std::map<std::string, ChemicalGroup>::const_iterator CTerminusIterator =
        mChemicalGroups.find("-COOH");
    if (CTerminusIterator == mChemicalGroups.end())
    {
        throw ChemicalBasisException("The default -COOH C-terminus not found.");
    }
    return CTerminusIterator->second;
}

double ChemicalBasis::secondSolventBindEnergy() const
{
    return mSecondSolventBindEnergy;
}

void ChemicalBasis::setSecondSolventBindEnergy(double newEnergy)
{
    mSecondSolventBindEnergy = newEnergy;
}

double ChemicalBasis::segmentLength() const
{
    return mSegmentLength;
}

void ChemicalBasis::setSegmentLength(double newSegmentLength)
{
    if (newSegmentLength <= 0.0)
    {
        throw ChemicalBasisException(
            "The new length of a segment is not positive.");
    }

    mSegmentLength = newSegmentLength;
}

int ChemicalBasis::kuhnLength() const
{
    return mKuhnLength;
}

void ChemicalBasis::setKuhnLength(int newKuhnLength)
{
    if (newKuhnLength <= 0)
    {
        throw ChemicalBasisException(
            "The new persistent length is not positive.");
    }
    mKuhnLength = newKuhnLength;
}

double ChemicalBasis::adsorbtionLayerWidth() const
{
    return mAdsorbtionLayerWidth;
}

void ChemicalBasis::setAdsorbtionLayerWidth(double newAdsorbtionLayerWidth)
{
    if (newAdsorbtionLayerWidth < 0.0)
    {
        throw ChemicalBasisException(
            "The new adsorbtion layer width is negative.");
    }
    mAdsorbtionLayerWidth = newAdsorbtionLayerWidth;
}

void ChemicalBasis::addChemicalGroup(ChemicalGroup newChemicalGroup)
{
    mChemicalGroups[newChemicalGroup.label()] = newChemicalGroup;
}

void ChemicalBasis::removeChemicalGroup(std::string label)
{
    if (mChemicalGroups.erase(label) == (std::size_t)0)
    {
        throw ChemicalBasisException(
            "The chemical group " + label + " is not found.");
    }
}

void ChemicalBasis::clearChemicalGroups()
{
    mChemicalGroups.clear();
}

void ChemicalBasis::setChemicalGroupBindEnergy(std::string label,
        double newBindEnergy)
{
    std::map<std::string,ChemicalGroup>::iterator it =
        mChemicalGroups.find(label);
    if (it == mChemicalGroups.end())
    {
        throw ChemicalBasisException(
            "The chemical group " + label + " is not found.");
    }
    it->second.setBindEnergy(newBindEnergy);
}

const ModelType ChemicalBasis::model() const
{
    return mModel;
}

void ChemicalBasis::setModel(ModelType newModel)
{
    mModel = newModel;
}

ChemicalBasis ChemicalBasis::setPredefinedChemicalBasis(
    PredefinedChemicalBasis predefinedChemicalBasisId)
{
    switch ( predefinedChemicalBasisId ) 
    {
        case RP_ACN_TFA_COIL_BOLTZMANN:
        {
            setModel(COIL_BOLTZMANN);
            clearChemicalGroups();
            addChemicalGroup(ChemicalGroup ("Alanine",
                                            "A",
                                            1.1425,
                                            71.0788,
                                            71.03711));
            addChemicalGroup(ChemicalGroup ("Cysteine",
                                            "C",
                                            1.2955,
                                            103.1388,
                                            103.00919));
            addChemicalGroup(ChemicalGroup ("Carboxyamidomethylated cysteine",
                                            "camC",
                                            0.77,
                                            160.1901,
                                            160.03065));
            addChemicalGroup(ChemicalGroup ("Aspartic_acid",
                                            "D",
                                            0.7805,
                                            115.0886,
                                            115.02694));
            addChemicalGroup(ChemicalGroup ("Glutamic_acid",
                                            "E",
                                            0.9835,
                                            129.1155,
                                            129.04259));
            addChemicalGroup(ChemicalGroup ("Phenylalanine",
                                            "F",
                                            2.3185,
                                            147.1766,
                                            147.06841));
            addChemicalGroup(ChemicalGroup ("Glycine",
                                            "G",
                                            0.6555,
                                            57.0519,
                                            57.02146));
            addChemicalGroup(ChemicalGroup ("Histidine",
                                            "H",
                                            0.3855,
                                            137.1411,
                                            137.05891));
            addChemicalGroup(ChemicalGroup ("Isoleucine",
                                            "I",
                                            2.1555,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Lysine",
                                            "K",
                                            0.2655,
                                            128.1741,
                                            128.09496));
            addChemicalGroup(ChemicalGroup ("Leucine",
                                            "L",
                                            2.2975,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Methionine",
                                            "M",
                                            1.8215,
                                            131.1926,
                                            131.04049));
            addChemicalGroup(ChemicalGroup ("Oxidated methionine",
                                            "oxM",
                                            1.8215,
                                            131.1926 + 15.9994,
                                            131.04049 + 15.994915));
            addChemicalGroup(ChemicalGroup ("Asparagine",
                                            "N",
                                            0.6135,
                                            114.1038,
                                            114.04293));
            addChemicalGroup(ChemicalGroup ("Proline",
                                            "P",
                                            1.1425,
                                            97.1167,
                                            97.05276));
            addChemicalGroup(ChemicalGroup ("Glutamine",
                                            "Q",
                                            0.7455,
                                            128.1307,
                                            128.05858));
            addChemicalGroup(ChemicalGroup ("Arginine",
                                            "R",
                                            0.5155,
                                            156.1875,
                                            156.10111));
            addChemicalGroup(ChemicalGroup ("Serine",
                                            "S",
                                            0.6975,
                                            87.0782,
                                            87.03203));
            addChemicalGroup(ChemicalGroup ("Phosphorylated serine",
                                            "pS",
                                            0.45,
                                            167.0581,
                                            166.99836));
            addChemicalGroup(ChemicalGroup ("Threonine",
                                            "T",
                                            0.8755,
                                            101.1051,
                                            101.04768));
            addChemicalGroup(ChemicalGroup ("Phosphorylated threonine",
                                            "pT",
                                            0.74,
                                            181.085,
                                            181.01401));
            addChemicalGroup(ChemicalGroup ("Valine",
                                            "V",
                                            1.7505,
                                            99.1326,
                                            99.06841));
            addChemicalGroup(ChemicalGroup ("Tryptophan",
                                            "W",
                                            2.4355,
                                            186.2132,
                                            186.07931));
            addChemicalGroup(ChemicalGroup ("Tyrosine",
                                            "Y",
                                            1.6855,
                                            163.176,
                                            163.06333));
            addChemicalGroup(ChemicalGroup ("Phosphorylated tyrosine",
                                            "pY",
                                            1.32,
                                            243.1559,
                                            243.02966));
            addChemicalGroup(ChemicalGroup("N-terminal hydrogen",
                                           "H-",
                                           -1.69,
                                           1.0079,
                                           1.00782));
            addChemicalGroup(ChemicalGroup("N-terminal acetyl",
                                           "Ac-",
                                           0.0,
                                           43.0452,
                                           43.01839));
            addChemicalGroup(ChemicalGroup("C-terminal carboxyl",
                                           "-COOH",
                                           -0.03,
                                           17.0073,
                                           17.00274));
            addChemicalGroup(ChemicalGroup("C-terminal amide",
                                           "-NH2",
                                           0.0,
                                           16.0226,
                                           16.01872));

            setSecondSolventBindEnergy(2.3979);
            setSegmentLength(10.0);
            setAdsorbtionLayerWidth(15.0);
            setKuhnLength(1);
            break;
        }

        case RP_ACN_FA_ROD_BOLTZMANN:
        {
            setModel(ROD_BOLTZMANN);
            clearChemicalGroups();
            addChemicalGroup(ChemicalGroup ("Alanine",
                                            "A",
                                            0.81,
                                            71.0788,
                                            71.03711));
            addChemicalGroup(ChemicalGroup ("Cysteine",
                                            "C",
                                            0.90,
                                            103.1388,
                                            103.00919));
            addChemicalGroup(ChemicalGroup ("Carboxyamidomethylated cysteine",
                                            "camC",
                                            0.24,
                                            160.1901,
                                            160.03065));
            addChemicalGroup(ChemicalGroup ("Aspartic_acid",
                                            "D",
                                            0.60,
                                            115.0886,
                                            115.02694));
            addChemicalGroup(ChemicalGroup ("Glutamic_acid",
                                            "E",
                                            0.62,
                                            129.1155,
                                            129.04259));
            addChemicalGroup(ChemicalGroup ("Phenylalanine",
                                            "F",
                                            2.65,
                                            147.1766,
                                            147.06841));
            addChemicalGroup(ChemicalGroup ("Glycine",
                                            "G",
                                            0.47,
                                            57.0519,
                                            57.02146));
            addChemicalGroup(ChemicalGroup ("Histidine",
                                            "H",
                                            -0.77,
                                            137.1411,
                                            137.05891));
            addChemicalGroup(ChemicalGroup ("Isoleucine",
                                            "I",
                                            2.20,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Lysine",
                                            "K",
                                            -0.64,
                                            128.1741,
                                            128.09496));
            addChemicalGroup(ChemicalGroup ("Leucine",
                                            "L",
                                            2.38,
                                            113.1594,
                                            113.08406));
            addChemicalGroup(ChemicalGroup ("Methionine",
                                            "M",
                                            1.73,
                                            131.1926,
                                            131.04049));
            addChemicalGroup(ChemicalGroup ("Oxidated methionine",
                                            "oxM",
                                            1.73,
                                            131.1926 + 15.9994,
                                            131.04049 + 15.994915));
            addChemicalGroup(ChemicalGroup ("Asparagine",
                                            "N",
                                            0.28,
                                            114.1038,
                                            114.04293));
            addChemicalGroup(ChemicalGroup ("Proline",
                                            "P",
                                            0.62,
                                            97.1167,
                                            97.05276));
            addChemicalGroup(ChemicalGroup ("Glutamine",
                                            "Q",
                                            0.47,
                                            128.1307,
                                            128.05858));
            addChemicalGroup(ChemicalGroup ("Arginine",
                                            "R",
                                            -0.59,
                                            156.1875,
                                            156.10111));
            addChemicalGroup(ChemicalGroup ("Serine",
                                            "S",
                                            0.50,
                                            87.0782,
                                            87.03203));
            addChemicalGroup(ChemicalGroup ("Phosphorylated serine",
                                            "pS",
                                            0.55,
                                            167.0581,
                                            166.99836));
            addChemicalGroup(ChemicalGroup ("Threonine",
                                            "T",
                                            0.69,
                                            101.1051,
                                            101.04768));
            addChemicalGroup(ChemicalGroup ("Phosphorylated threonine",
                                            "pT",
                                            0.59,
                                            181.085,
                                            181.01401));
            addChemicalGroup(ChemicalGroup ("Valine",
                                            "V",
                                            1.49,
                                            99.1326,
                                            99.06841));
            addChemicalGroup(ChemicalGroup ("Tryptophan",
                                            "W",
                                            2.87,
                                            186.2132,
                                            186.07931));
            addChemicalGroup(ChemicalGroup ("Tyrosine",
                                            "Y",
                                            1.40,
                                            163.176,
                                            163.06333));
            addChemicalGroup(ChemicalGroup ("Phosphorylated tyrosine",
                                            "pY",
                                            0.92,
                                            243.1559,
                                            243.02966));
            addChemicalGroup(ChemicalGroup("N-terminal hydrogen",
                                           "H-",
                                           -2.39,
                                           1.0079,
                                           1.00782));
            addChemicalGroup(ChemicalGroup("N-terminal acetyl",
                                           "Ac-",
                                           -0.05,
                                           43.0452,
                                           43.01839));
            addChemicalGroup(ChemicalGroup("C-terminal carboxyl",
                                           "-COOH",
                                           1.02,
                                           17.0073,
                                           17.00274));
            addChemicalGroup(ChemicalGroup("C-terminal amide",
                                           "-NH2",
                                           0.55,
                                           16.0226,
                                           16.01872));

            setSecondSolventBindEnergy(2.4);
            setSegmentLength(4.0);
            setAdsorbtionLayerWidth(16.0);
            setKuhnLength(1);
            break;
        }
    }

    return (*this);
}
}
