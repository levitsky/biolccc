#include "chemicalbasis.h"

namespace BioLCCC {

//TODO: make addAminoacid method!
ChemicalBasis::ChemicalBasis() {
    
    //adding standard aminoacids, masses set to zero.
    addAminoacid(Aminoacid ("Alanine",
                                "A",
                                1.1425,
                                71.0788,
                                71.03711));
    addAminoacid(Aminoacid ("Cysteine",
                                "C",
                                1.2955,
                                103.1388,
                                103.00919));
    addAminoacid(Aminoacid ("Carboxyamidomethylated cysteine",
                                "camC",
                                0.77,
                                160.1901,
                                160.03065));
    addAminoacid(Aminoacid ("Aspartic_acid",
                                "D",
                                0.7805,
                                115.0886,
                                115.02694));
    addAminoacid(Aminoacid ("Glutamic_acid",
                                "E",
                                0.9835,
                                129.1155,
                                129.04259));
    addAminoacid(Aminoacid ("Phenylalanine",
                                "F",
                                2.3185,
                                147.1766,
                                147.06841));
    addAminoacid(Aminoacid ("Glycine",
                                "G",
                                0.6555,
                                57.0519,
                                57.02146));
    addAminoacid(Aminoacid ("Histidine",
                                "H",
                                0.3855,
                                137.1411,
                                137.05891));
    addAminoacid(Aminoacid ("Isoleucine",
                                "I",
                                2.1555,
                                113.1594,
                                113.08406));
    addAminoacid(Aminoacid ("Lysine",
                                "K",
                                0.2655,
                                128.1741,
                                128.09496));
    addAminoacid(Aminoacid ("Leucine",
                                "L",
                                2.2975,
                                113.1594,
                                113.08406));
    addAminoacid(Aminoacid ("Methionine",
                                "M",
                                1.8215,
                                131.1926,
                                131.04049));
    addAminoacid(Aminoacid ("Asparagine",
                                "N",
                                0.6135,
                                114.1038,
                                114.04293));
    addAminoacid(Aminoacid ("Proline",
                                "P",
                                1.1425,
                                97.1167,
                                97.05276));
    addAminoacid(Aminoacid ("Glutamine",
                                "Q",
                                0.7455,
                                128.1307,
                                128.05858));
    addAminoacid(Aminoacid ("Arginine",
                                "R",
                                0.5155,
                                156.1875,
                                156.10111));
    addAminoacid(Aminoacid ("Serine",
                                "S",
                                0.6975,
                                87.0782,
                                87.03203));
    addAminoacid(Aminoacid ("Phosphorylated serine",
                                "pS",
                                0.45,
                                167.0581,
                                166.99836));
    addAminoacid(Aminoacid ("Threonine",
                                "T",
                                0.8755,
                                101.1051,
                                101.04768));
    addAminoacid(Aminoacid ("Phosphorylated threonine",
                                "pT",
                                0.74,
                                181.085,
                                181.01401));
    addAminoacid(Aminoacid ("Valine",
                                "V",
                                1.7505,
                                99.1326,
                                99.06841));
    addAminoacid(Aminoacid ("Tryptophan",
                                "W",
                                2.4355,
                                186.2132,
                                186.07931));
    addAminoacid(Aminoacid ("Tyrosine",
                                "Y",
                                1.6855,
                                163.176,
                                163.06333));
    addAminoacid(Aminoacid ("Phosphorylated tyrosine",
                                "pY",
                                1.32,
                                243.1559,
                                243.02966));
                                    
    //adding standard terminal groups
    addNTerminus(Terminus("N-terminal hydrogen",
                            "H-",
                            -1.69,
                            1.0079,
                            1.00782));
    addNTerminus(Terminus("N-terminal acetyl",
                            "Ac-",
                            0,
                            43.0452,
                            43.01839));
    addCTerminus(Terminus("C-terminal carboxyl",
                            "-COOH",
                            -0.03,
                            17.0073,
                            17.00274));
    addCTerminus(Terminus("C-terminal amid",
                            "-NH2",
                            0,
                            16.0226,
                            16.01872));
                                    
    //setting standard second solvent bind energy 
    mSecondSolventBindEnergy = 2.3979;
}

const std::map<std::string,Aminoacid> & ChemicalBasis::aminoacids() const{
    return mAminoacids;
}

const std::map<std::string,Terminus> & ChemicalBasis::NTermini() const{
    return mNTermini;
}

const std::map<std::string,Terminus> & ChemicalBasis::CTermini() const{
    return mCTermini;
}

const Terminus & ChemicalBasis::defaultNTerminus() const {
    return mNTermini.find("H-")->second;
}

const Terminus & ChemicalBasis::defaultCTerminus() const{
    return mCTermini.find("-COOH")->second;
}

double ChemicalBasis::secondSolventBindEnergy() const{
    return mSecondSolventBindEnergy;
}

void ChemicalBasis::setSecondSolventBindEnergy(double newEnergy){
    mSecondSolventBindEnergy = newEnergy;
}

void ChemicalBasis::addAminoacid(Aminoacid newAminoacid){
    mAminoacids[newAminoacid.label()] = newAminoacid;
}

void ChemicalBasis::addNTerminus(Terminus newNTerminus){
    mNTermini[newNTerminus.label()] = newNTerminus;
}

void ChemicalBasis::addCTerminus(Terminus newCTerminus){
    mCTermini[newCTerminus.label()] = newCTerminus;
}

bool ChemicalBasis::removeAminoacid(std::string label){
    return (bool) mAminoacids.erase(label);
}

bool ChemicalBasis::removeNTerminus(std::string label){
    return (bool) mNTermini.erase(label);
}

bool ChemicalBasis::removeCTerminus(std::string label){
    return (bool) mCTermini.erase(label);
}

bool ChemicalBasis::setAminoacidBindEnergy(std::string label, double newBindEnergy){
    std::map<std::string,Aminoacid>::iterator it = mAminoacids.find(label);
    if (it == mAminoacids.end()) return false;
    it->second.setBindEnergy(newBindEnergy);
    return true;
}

bool ChemicalBasis::setNTerminusBindEnergy(std::string label, double newBindEnergy){
    std::map<std::string,Terminus>::iterator it = mNTermini.find(label);
    if (it == mNTermini.end()) return false;
    it->second.setBindEnergy(newBindEnergy);
    return true;
}

bool ChemicalBasis::setCTerminusBindEnergy(std::string label, double newBindEnergy){
    std::map<std::string,Terminus>::iterator it = mCTermini.find(label);
    if (it == mCTermini.end()) return false;
    it->second.setBindEnergy(newBindEnergy);
    return true;
}
}
