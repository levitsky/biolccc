#ifndef AMINOACID_H
#define AMINOACID_H

#include <string>

//! The Aminoacid class encapsulates all physical properties of an aminoacid.

namespace BioLCCC {

class Aminoacid {
    
    public:
        /*!
            Constructs an aminoacid with given Name, Label, bind energy and 
            masses.
        */
        Aminoacid(std::string name = "",
                std::string label = "",
                double bindEnergy = 0.0,
                double averageMass = 0.0,
                double monoisotopicMass = 0.0
        );
        
        /*!
            Returns the full name of the aminoacid.
            
            Examples:
            - Methionine
            - Phosphorylated threonine
        */
        std::string name() const;
        
        /*!
            Returns the brief code of the aminoacid, which is used in peptide 
            sequence notation.
            
            Examples:
            - M
            - pT
            
            The set of IUPAC conventions for peptide sequence notation could be
            found at http://www.chem.qmul.ac.uk/iupac/AminoAcid/.
        */
        std::string label() const;
        
        /*!
            Returns the average mass of the aminoacid residue.
            
            The average mass is measured for R-CH(NH-)-CO- structure WITHOUT 
            terminal H- and -OH (equals to monoisotopic mass of a whole 
            aminoacid molecule minus 18.01528).
        */
        double averageMass() const;
        
        /*!
            Returns the monoisotopic mass of the aminoacid residue.
        
            The monoisotopic mass is measured for R-CH(NH-)-CO- structure 
            WITHOUT terminal H- and -OH (equals to monoisotopic mass of 
            a whole aminoacid molecule minus 18.010565).
        */
        double monoisotopicMass() const;
        
        /*!
            Returns the bind energy of the aminoacid according to the BioLCCC 
            model. The zero is bind energy of water and the unit is kT.
        
            The current version of BioLCCC contains only one energy per 
            aminoacid: the energy of interaction between reversed phase C18 
            (octadecyl) and in-chain aminoacid residue at 293 K and pH = 2.0. 
            We belive that the next versions of the BioLCCC model will contain 
            energies for other chromatographic phases and conditions, so this 
            part of code will be rewrited.
        */
        double bindEnergy() const;
    
    private:
        std::string mName;
        std::string mLabel;
        double mBindEnergy;
        double mAverageMass;
        double mMonoisotopicMass;
};

}

#endif
