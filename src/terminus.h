#ifndef TERMINUS_H
#define TERMINUS_H

#include <string>

//! The Terminus class encapsulates all physical properties of a peptide 
//! terminal groups.


namespace BioLCCC {
class Terminus {
    
    public:
        
        /*!
            Constructs a peptide terminal group with given name, label, 
            bind energy and masses.
        */
        Terminus(std::string name = "",
                std::string label = "",
                double bindEnergy = 0.0,
                double averageMass = 0.0,
                double monoisotopicMass = 0.0
        );
        
        /*!
            Returns a full name of the peptide terminal group.
            
            Examples:
            - N-terminal acetyl
            - C-terminal carboxyl
        */
        std::string name() const;
        
        /*!
            Returns the brief code of the peptide terminal group, which is used
            in peptide sequence notation.
            
            Examples:
            - Ac-
            - -COOH
            
            The set of IUPAC conventions for peptide sequence notation could be
            found at http://www.chem.qmul.ac.uk/iupac/AminoAcid/.
        */
        std::string label() const;
        
        /*!
            Returns the average mass of the terminal group.
        */
        double averageMass() const;
        
        /*!
            Returns the monoisotopic mass of the terminal group.
        */
        double monoisotopicMass() const;
        
        /*!
            Returns the modification of bind energy for the aminoacid with 
            this terminal group according to the BioLCCC model. 
            The zero is bind energy of water and the unit is kT.
        
            The current version of BioLCCC contains only one energy per
            aminoacid: the energy of interaction between reversed phase C18 
            (octadecyl) and an in-chain aminoacid residue at 293 K and pH = 2.0.
            We hope that next versions of the BioLCCC model will contain 
            energies for other chromatographic phases and conditions, so this 
            part of code will be rewrited.
        */
        double bindEnergy() const;
    
    private:
        std::string mName;
        std::string mLabel;
        double mAverageMass;
        double mMonoisotopicMass;
        double mBindEnergy;
};

}

#endif
