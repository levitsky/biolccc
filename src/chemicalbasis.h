#ifndef CHEMICALBASIS_H
#define CHEMICALBASIS_H

#include <map>
#include "aminoacid.h"
#include "terminus.h"

namespace BioLCCC {

//! The ChemicalBasis class manages aminoacids and terminal groups.
/*!
    <b>Discussing the internal structure of a ChemicalBasis class.</b>
    There are two possible ways to store aminoacids here: a map and a vector.
    - A map is simple in realization. We could write a public function
    getAminoacid(std::string) and use it to parse peptide sequence. But in this
    case we have a problem with two-and-more-letters labels. Actually, we need 
    to have the full list of aminoacids to be sure, whether peptide sequence 
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

class ChemicalBasis {
    
    public:
        
        /*!
            Constructs a new ChemicalBasis object with the standard BioLCCC list
            of aminoacids, termini and the standard energy of second solvent.
        */
        ChemicalBasis();
        
        /*!
            Returns the map of aminoacids.
        */
        const std::map<std::string,Aminoacid> & aminoacids () const;
        
        /*!
            Returns the map of N-Terminal groups.
        */
        const std::map<std::string,Terminus> & NTermini() const;
        
        /*!
            Returns an iterator to the end of the C-Terminal group  map.
        */
        const std::map<std::string,Terminus> & CTermini() const;
        
        /*!
            Returns an iterator to the default N-terminal group.
        */
        const Terminus & defaultNTerminus() const;

        /*!
            Returns an iterator to the default N-terminal group.
        */
        const Terminus & defaultCTerminus() const;
        
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
            the number of aminoacids between the joints of a polymer.
        */
        bool setPersistentLength(int newPersistentLength);

        /*!
            Returns the persistent length of a biopolymer. Persistent length 
            equals the number of aminoacids between the joints of a polymer.
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
            Adds a new aminoacid.
        */
        void addAminoacid(Aminoacid newAminoacid);

        /*!
            Adds a new N-terminus.
        */
        void addNTerminus(Terminus newNTerminus);

        /*!
            Adds a new C-terminus.
        */
        void addCTerminus(Terminus newCTerminus);

        /*!
            Removes the aminoacid with the given label; returns 'true' on 
            success.
        */
        bool removeAminoacid (std::string label);

        /*!
            Removes the N-terminus with the given label; returns 'true' on 
            success.
        */
        bool removeNTerminus (std::string label);

        /*!
            Removes the C-terminus with the given label; returns 'true' on 
            success.
        */
        bool removeCTerminus (std::string label);

        /*!
            Sets the value of binding energy for the aminoacid with the given 
            label;
            returns 'true' on success or 'false' if the aminoacid is not found.
        */
        bool setAminoacidBindEnergy (std::string label, double newBindEnergy);

        /*!
            Sets the value of binding energy for the N-terminus with the given 
            label;
            returns 'true' on success or 'false' if the aminoacid is not found.
        */
        bool setNTerminusBindEnergy (std::string label, double newBindEnergy);

        /*!
            Sets the value of binding energy for the C-terminus with the given 
            label;
            returns 'true' on success or 'false' if the aminoacid is not found.
        */
        bool setCTerminusBindEnergy (std::string label, double newBindEnergy);

        /*!
            Sets the type of BioLCCC model (e.g. CoilBoltzmann, CoilSnyder,
            RodBoltzmann and so on).
        */
        bool setModel(std::string newModel);

        /*!
            Returns the type of BioLCCC model which is used in calculations with
            this ChemicalBasis.
        */
        const std::string model() const;

    private:
        std::map<std::string,Aminoacid> mAminoacids;
        std::map<std::string,Terminus>  mNTermini;
        std::map<std::string,Terminus>  mCTermini;
        double mSecondSolventBindEnergy;
        double mSegmentLength;
        int mPeristentLength;
        double mAdsorbtionLayerWidth;
        std::string mModel;
};

}

#endif
