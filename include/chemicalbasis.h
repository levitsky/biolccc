#ifndef CHEMICALBASIS_H
#define CHEMICALBASIS_H

#include <map>
#include "chemicalgroup.h"

namespace BioLCCC {

//! The ChemicalBasis class manages amino acids and terminal groups.
/*!
    <b>Discussing the internal structure of a ChemicalBasis class.</b>
    There are two possible ways to store amino acids here: a map and a vector.
    - A map is simple in realization. We could write a public function
    getAminoacid(std::string) and use it to parse peptide sequence. But in this
    case we have a problem with two-and-more-letters labels. Actually, we need 
    to have the full list of amino acids to be sure, whether peptide sequence 
    contains an unidentifiable label, or we just need to check more letters 
    from the sequence.
    - A vector needs to implement two functions, beginIterator and endIterator
    to incapsulate all data. But in that case we could easily check the whole 
    list at one time.
    BUT! I've found that map could be used as a vector too, so we implement 
    a map.
    <b>Default terminal group</b>
    Be careful, the default terminal groups should be the first in the 
    corresponding vector!
*/

enum ModelType {
    COIL_BOLTZMANN,
    COIL_BOLTZMANN_DOUBLE_LAYER,
    COIL_SNYDER,
    ROD_BOLTZMANN
};

class ChemicalBasis {
    public:
        
        /*!
            Constructs a new ChemicalBasis object with the standard BioLCCC list
            of amino acids, termini and the standard energy of second solvent.
        */
        ChemicalBasis();
        
        /*!
            Returns the map of all chemical groups.
        */
        const std::map<std::string, ChemicalGroup> & chemicalGroups() const;
        
        /*!
            Returns an iterator to the default N-terminal group.
        */
        const ChemicalGroup & defaultNTerminus() const;

        /*!
            Returns an iterator to the default N-terminal group.
        */
        const ChemicalGroup & defaultCTerminus() const;
        
        /*!
            Returns the bind energy of the second solvent. The zero is the bind 
            energy of water and the unit is kT.
        */
        double secondSolventBindEnergy() const;
        
        /*!
            Sets a new value of the bind energy for the second solvent. The zero
            is bind energy of water and the unit is kT.
        */
        void setSecondSolventBindEnergy(double newEnergy);

        /*!
            Return the length between two peptide bonds, angstroms.
        */
        double segmentLength() const;

        /*!
            Sets the length between two peptide bonds, angstroms.
        */
        bool setSegmentLength(double newSegmentLength);

        /*!
            Sets the persistent length of a biopolymer. Persistent length equals
            the number of amino acids between the joints of a polymer.
        */
        bool setPersistentLength(int newPersistentLength);

        /*!
            Returns the persistent length of a biopolymer. Persistent length 
            equals the number of amino acids between the joints of a polymer.
        */
        int persistentLength() const;

        /*!
            Return the width of a solid phase adsorbtion layer.
        */
        double adsorbtionLayerWidth() const;

        /*!
            Sets a new value of the width of a solid phase adsorbtion layer.
        */
        bool setAdsorbtionLayerWidth(double newAdsorbtionLayerWidth);

        /*!
            Adds a new chemical group.
        */
        void addChemicalGroup(ChemicalGroup newChemicalGroup);

        /*!
            Removes the chemical group with the given label; returns 'true' on 
            success.
        */
        bool removeChemicalGroup(std::string label);

        /*!
            Removes all chemical groups in a basis.
        */
        void clearChemicalGroups();

        /*!
            Sets the value of binding energy for the chemical group with the 
            given label;
            Returns 'true' on success or 'false' if the chemical group is not 
            found.
        */
        bool setChemicalGroupBindEnergy(std::string label,double newBindEnergy);

        /*!
            Sets the type of BioLCCC model (e.g. CoilBoltzmann, CoilSnyder,
            RodBoltzmann and so on).
        */
        bool setModel(ModelType newModel);

        /*!
            Returns the type of BioLCCC model which is used in calculations with
            this ChemicalBasis.
        */
        const ModelType model() const;

    private:
        std::map<std::string,ChemicalGroup> mChemicalGroups;
        double mSecondSolventBindEnergy;
        double mSegmentLength;
        int mPeristentLength;
        double mAdsorbtionLayerWidth;
        ModelType mModel;
};

}

#endif
