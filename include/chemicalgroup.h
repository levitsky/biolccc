#ifndef AMINOACID_H
#define AMINOACID_H

#include <string>

/*!
    The ChemicalGroup class encapsulates all physical properties of a group of
    atoms. This group may be an amino acid residue or a terminal group,
    depending on its label. 
*/

namespace BioLCCC {

class ChemicalGroup {
    
    public:
        /*!
            Constructs a chemical group with given Name, Label, bind energy and 
            masses.
        */
        ChemicalGroup(std::string name = "",
                std::string label = "",
                double bindEnergy = 0.0,
                double averageMass = 0.0,
                double monoisotopicMass = 0.0
        );
        
        /*!
            Returns the full name of the chemical group.
            
            Examples:
            - Methionine
            - Phosphorylated threonine
            - C-Terminal amidation
        */
        std::string name() const;
        
        /*!
            Returns the brief code of the group, which is used in peptide 
            sequence notation. For an amino acid it is one-letter code in the
            case of the standard amino acids, or extended code for the modified
            aminoacids. For a terminal group this code have to start or end with
            a dash.
            
            Examples:
            - M
            - pT
            - Ac-
            
            The set of IUPAC conventions for peptide sequence notation could be
            found at http://www.chem.qmul.ac.uk/iupac/AminoAcid/.
        */
        std::string label() const;
        
        /*!
            Returns the average mass of the chemical group.
            
            The average mass of an amino acid is measured for R-CH(NH-)-CO- 
            structure WITHOUT terminal H- and -OH (equals to the average mass 
            of a whole aminoacid molecule minus 18.01528).
        */
        double averageMass() const;
        
        /*!
            Returns the monoisotopic mass of the chemical group.
        
            The monoisotopic mass of an amino acid is measured for R-CH(NH-)-CO-
            structure WITHOUT terminal H- and -OH (equals to the monoisotopic 
            mass of a whole amino acid molecule minus 18.010565).
        */
        double monoisotopicMass() const;
        
        /*!
            Returns the bind energy of the chemical group according to the 
            BioLCCC model. The zero is bind energy of water and the unit is kT.
       
            The bind energy of a terminal group is added to the binding group of
            the corresponding terminal amino acid.
        */
        double bindEnergy() const;

        /*!
            Shows whether the group is N-Terminal.
        */
        bool isNTerminal() const;

        /*!
            Shows whether the group is C-Terminal.
        */
        bool isCTerminal() const;

        /*!
            Shows whether the group is an amino acid.
        */
        bool isAminoAcid() const;

        /*!
            Set the binding energy value.
        */
        void setBindEnergy(double newBindEnergy);

    private:
        std::string mName;
        std::string mLabel;
        double mBindEnergy;
        double mAverageMass;
        double mMonoisotopicMass;
};

}

#endif
