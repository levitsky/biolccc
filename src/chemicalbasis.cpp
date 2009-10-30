#include "chemicalbasis.h"

namespace BioLCCC {

ChemicalBasis::ChemicalBasis() {
    
    //adding standard aminoacids, masses set to zero.
    mAminoacids["A"] = Aminoacid ("Alanine",
                                "A",
                                1.1425,
                                71.0788,
                                71.03711);
    mAminoacids["C"] = Aminoacid ("Cysteine",
                                "C",
                                1.2955,
                                103.1388,
                                103.00919);
    mAminoacids["D"] = Aminoacid ("Aspartic_acid",
                                "D",
                                0.7805,
                                115.0886,
                                115.02694);
    mAminoacids["E"] = Aminoacid ("Glutamic_acid",
                                "E",
                                0.9835,
                                129.1155,
                                129.04259);
    mAminoacids["F"] = Aminoacid ("Phenylalanine",
                                "F",
                                2.3185,
                                147.1766,
                                147.06841);
    mAminoacids["G"] = Aminoacid ("Glycine",
                                "G",
                                0.6555,
                                57.0519,
                                57.02146);
    mAminoacids["H"] = Aminoacid ("Histidine",
                                "H",
                                0.3855,
                                137.1411,
                                137.05891);
    mAminoacids["I"] = Aminoacid ("Isoleucine",
                                "I",
                                2.1555,
                                113.1594,
                                113.08406);
    mAminoacids["K"] = Aminoacid ("Lysine",
                                "K",
                                0.2655,
                                128.1741,
                                128.09496);
    mAminoacids["L"] = Aminoacid ("Leucine",
                                "L",
                                2.2975,
                                113.1594,
                                113.08406);
    mAminoacids["M"] = Aminoacid ("Methionine",
                                "M",
                                1.8215,
                                131.1926,
                                131.04049);
    mAminoacids["N"] = Aminoacid ("Asparagine",
                                "N",
                                0.6135,
                                114.1038,
                                114.04293);
    mAminoacids["P"] = Aminoacid ("Proline",
                                "P",
                                1.1425,
                                97.1167,
                                97.05276);
    mAminoacids["Q"] = Aminoacid ("Glutamine",
                                "Q",
                                0.7455,
                                128.1307,
                                128.05858);
    mAminoacids["R"] = Aminoacid ("Arginine",
                                "R",
                                0.5155,
                                156.1875,
                                156.10111);
    mAminoacids["S"] = Aminoacid ("Serine",
                                "S",
                                0.6975,
                                87.0782,
                                87.03203);
    mAminoacids["T"] = Aminoacid ("Threonine",
                                "T",
                                0.8755,
                                101.1051,
                                101.04768);
    mAminoacids["V"] = Aminoacid ("Valine",
                                "V",
                                1.7505,
                                99.1326,
                                99.06841);
    mAminoacids["W"] = Aminoacid ("Tryptophan",
                                "W",
                                2.4355,
                                186.2132,
                                186.07931);
    mAminoacids["Y"] = Aminoacid ("Tyrosine",
                                "Y",
                                1.6855,
                                163.176,
                                163.06333);
                                    
    //adding standard terminal groups
    mNTermini["H-"] = Terminus ("N-terminal hydrogen",
                            "H-",
                            -1.69,
                            1.0079,
                            1.00782);
    mNTermini["Ac-"] = Terminus ("N-terminal acetyl",
                            "Ac-",
                            0,
                            43.0452,
                            43.01839);
    mCTermini["-COOH"] = Terminus ("C-terminal carboxyl",
                            "-COOH",
                            -0.03,
                            17.0073,
                            17.00274);
    mCTermini["-NH2"] = Terminus ("C-terminal amid",
                            "-NH2",
                            0,
                            16.0226,
                            16.01872);
                                    
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

void ChemicalBasis::setSecondSolventBindEnergy(double newEnergy) {
    mSecondSolventBindEnergy = newEnergy;
}

}
