#ifndef CHEMICALBASIS_H
#define CHEMICALBASIS_H

#include <map>
#include "biolcccexception.h"
#include "chemicalgroup.h"

namespace BioLCCC
{

//! This exception is raised when something goes wrong with a ChemicalBasis.
class ChemicalBasisException : public BioLCCCException
{
public:
    //! Constructs an instance of ChemicalBasis with a given message.
    ChemicalBasisException(std::string message);
};

//! A set of assumptions of the BioLCCC model.
/*!
    Different models come from different sets of initial assumption. For the
    better explanation see the theory of the BioLCCC model.
 */
enum ModelType
{
    COIL_BOLTZMANN, /*!< The standard BioLCCC model with the assumption 
        of an absolute flexibility of a protein molecule. */
    ROD_BOLTZMANN, /*!< The BioLCCC model with the assumption of absolute 
        rigidity of a peptide molecule. It works better for short molecules. */
    COIL_BOLTZMANN_DOUBLE_LAYER, /*!< EXPERIMENTAL. Modification of the 
        standard model in which the adsorption occurs in a volume, 
        i.e. in two near-wall layers. */
    COIL_SNYDER /*! EXPERIMENTAL. Modification of the standard model in which
        the adsorption is described by the linear Snyder's theory. */
};

//! An instance of ChemicalBasis contains a set of BioLCCC constants.
/*!
    An instance of ChemicalBasis manages all the physicochemical constants,
    which are used in the calculations. Currently, it contains:
        - The list of amino acids and peptide terminal groups.
        - The terminal groups which are set by default (cannon be changed).
        - The Energy of binding between a solvent and the surface of a solid 
          phase.
        - The type of BioLCCC model being used in calculations. 
        - Peptide geometry: the length of amino acid and the Kuhn length.
        - The width of the adsorbing layer.
       
    Note, that the set of constants is highly interconnected. Usually the change
    in one constant, like the width of adsorbing layer or type of BioLCCC model,
    would only deteriorate the quality of RT prediction.
 */
class ChemicalBasis
{
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
    void setSegmentLength(double newSegmentLength);

    /*!
        Sets the persistent length of a biopolymer. Persistent length equals
        the number of amino acids between the joints of a polymer.
    */
    void setKuhnLength(int newKuhnLength);

    /*!
        Returns the persistent length of a biopolymer. Persistent length
        equals the number of amino acids between the joints of a polymer.
    */
    int kuhnLength() const;

    /*!
        Return the width of a solid phase adsorbtion layer.
    */
    double adsorbtionLayerWidth() const;

    /*!
        Sets a new value of the width of a solid phase adsorbtion layer.
    */
    void setAdsorbtionLayerWidth(double newAdsorbtionLayerWidth);

    /*!
        Adds a new chemical group.
    */
    void addChemicalGroup(ChemicalGroup newChemicalGroup);

    /*!
        Removes the chemical group with the given label;
        Throws ChemicalBasisException if the chemical group is not found.
    */
    void removeChemicalGroup(std::string label);

    /*!
        Removes all chemical groups in a basis.
    */
    void clearChemicalGroups();

    /*!
        Sets the value of binding energy for the chemical group with the
        given label;
        Throws ChemicalBasisException if the chemical group is not found.
    */
    void setChemicalGroupBindEnergy(std::string label,double newBindEnergy);

    /*!
        Sets the type of BioLCCC model (e.g. CoilBoltzmann, CoilSnydel).
    */
    void setModel(ModelType newModel);

    /*!
        Returns the type of BioLCCC model which is used in calculations with
        this ChemicalBasis.
    */
    const ModelType model() const;

private:
    std::map<std::string,ChemicalGroup> mChemicalGroups;
    double mSecondSolventBindEnergy;
    double mSegmentLength;
    int mKuhnLength;
    double mAdsorbtionLayerWidth;
    ModelType mModel;
};

}

#endif
