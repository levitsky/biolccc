#ifndef CHEMICALBASIS_H
#define CHEMICALBASIS_H

#include <map>
#include "aminoacid.h"
#include "terminus.h"

//! The ChemicalBasis class manages aminoacids and terminal groups.
/*!
    <b>Discussing the internal structure of a ChemicalBasis class.</b>
    There are two possible ways to store aminoacids here: a map and a vector.
    - A map is simple in realization. We could write a public function getAminoacid(std::string) and use it to parse peptide sequence. But in this case we have a problem with two-and-more-letters labels. Actually, we need to have a full list of aminoacids to be sure, whether peptide sequence contains an unidentifiable label, or we just need to check more letters from a sequence.
    - A vector needs to implement two functions, beginIterator and endIterator to incapsulate all data. But in that case we could easily check the whole list at one time.
    BUT! I've found that map could be used as a vector too, so we implement a map.
    <b>Default terminal group</b>
    Be careful, the default terminal groups should be the first in corresponding vector!
*/

class ChemicalBasis {
    
    public:
        
        /*!
            Constructs a new ChemicalBasis object with a standard BioLCCC list of aminoacids, termini and standard energy of a second solvent.
        */
        ChemicalBasis();
        
        /*!
            Returns an iterator to the first aminoacid in the map.
        */
        std::map<std::string,Aminoacid>::const_iterator beginAminoacid() const;
        
        /*!
            Returns an iterator to the end of the aminoacid map.
        */
        std::map<std::string,Aminoacid>::const_iterator endAminoacid() const;
        
        /*!
            Returns an iterator to the first N-Terminal group in the map.
        */
        std::map<std::string,Terminus>::const_iterator beginNTerminus() const;
        
        /*!
            Returns an iterator to the end of the N-Terminal groups' map.
        */
        std::map<std::string,Terminus>::const_iterator endNTerminus() const;
        
        /*!
            Returns an iterator to the first C-Terminal group in the map.
        */
        std::map<std::string,Terminus>::const_iterator beginCTerminus() const;
        
        /*!
            Returns an iterator to the end of the C-Terminal group  map.
        */
        std::map<std::string,Terminus>::const_iterator endCTerminus() const;
        
        /*!
            Returns an iterator to default N-terminal group.
        */
        std::map<std::string,Terminus>::const_iterator defaultNTerminus() const;
        
        /*!
            Returns an iterator to default C-terminal group.
        */
        std::map<std::string,Terminus>::const_iterator defaultCTerminus() const;
        
        /*!
            Returns a bind energy of a second solvent. The zero is bind energy of water and the unit is kT.
        */
        double secondSolventBindEnergy() const;
        
        /*!
            Sets a new value of a bind energy for a second solvent. The zero is bind energy of water and the unit is kT.
        */
        void setSecondSolventBindEnergy(double newEnergy);
        
    private:
        std::map<std::string,Aminoacid> mAminoacids;
        std::map<std::string,Terminus>  mNTermini;
        std::map<std::string,Terminus>  mCTermini;
        double mSecondSolventBindEnergy;
};

#endif
