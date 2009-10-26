#ifndef AMINOACID_H
#define AMINOACID_H

#include <string>

//! The Aminoacid class encapsulates all physical properties of an aminoacid.


class Aminoacid {
    
    public:
        
        /*!
            Constructs an aminoacid with given Name, Label, bind energy and masses.
        */
        Aminoacid(std::string Name = "",
                std::string Label = "",
                double bindEnergy = 0.0,
                double averageMass = 0.0,
                double monoisotopicMass = 0.0
        );
        
        /*!
            Returns a full name of the aminoacid.
            
            Examples:
            - Methionine
            - Phosphorylated threonine
        */
        std::string name() const;
        
        /*!
            Returns a brief code of the aminoacid, which is used in peptide sequence notation.
            
            Examples:
            - M
            - pT
            
            A set of IUPAC conventions for a peptide sequence notation could be found at http://www.chem.qmul.ac.uk/iupac/AminoAcid/.
        */
        std::string label() const;
        
        /*!
            Returns an average mass of an aminoacid residue.
            
            The average mass is measured for R-CH(NH-)-CO- structure WITHOUT terminal H- and -OH (equals to monoisotopic mass of a whole aminoacid molecule minus 18.01528).
        */
        double averageMass() const;
        
        /*!
            Returns a monoisotopic mass of an aminoacid residue.
        
            The monoisotopic mass is measured for R-CH(NH-)-CO- structure WITHOUT terminal H- and -OH (equals to monoisotopic mass of a whole aminoacid molecule minus 18.010565).
        */
        double monoisotopicMass() const;
        
        /*!
            Returns a bind energy of an aminoacid according to the BioLCCC model. The zero is bind energy of water and the unit is kT.
        
            The current version of BioLCCC contains only one energy per aminoacid: an energy of interaction between reversed phase C18 (octadecyl) and in-chain aminoacid residue at 293 K and pH = 2.0. We hope that next versions of the BioLCCC model will contain energies for other chromatographic phases and conditions, so this part of code will be rewrited.
        */
        double bindEnergy() const;
    
    private:
        std::string mName;
        std::string mLabel;
        double mBindEnergy;
        double mAverageMass;
        double mMonoisotopicMass;
};

#endif