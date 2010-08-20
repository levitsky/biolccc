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
    //! Constructs an instance of ChemicalBasisException with the given message.
    ChemicalBasisException(std::string message);
};

//! A set of assumptions of the BioLCCC model.
/*!
    Different models come from different sets of initial assumption. For the
    better explanation see the theory of the BioLCCC model.
 */
enum ModelType
{
    /*! The standard BioLCCC model with the assumption of absolute 
        flexibility of a protein molecule. */
    COIL_BOLTZMANN,
    /*! The BioLCCC model with the assumption of absolute rigidity of 
        a protein molecule. It works better for short molecules. */
    ROD_BOLTZMANN,
    /*! EXPERIMENTAL. Modification of the standard model in which adsorption 
        occurs in a volume, i.e. in two near-wall layers. */
    COIL_BOLTZMANN_DOUBLE_LAYER, 
    /*! EXPERIMENTAL. Modification of the standard model in which adsorption 
        is described by the linear Snyder's theory. */
    COIL_SNYDER 
};

//! This enum describes the predefined sets of physicochemical constants.
/*!
    The BioLCCC library contains several predifined sets of physicochemical
    constants. Please note that usually changing only one parameter in a whole
    set of constants deteriorate the quality of prediction.
 */
enum PredefinedChemicalBasis
{
    //! Reversed phase, ACN, trifluoracetic acid, COIL_BOLTZMANN model.
    /*! A ChemicalBasis calibrated for reversed phase, ACN as a second solvent,
        0.1% TFA and COIL_BOLTZMANN type of BioLCCC model. The data was 
        obtained in Guo et al, Journal of Chromatography, 359 (1986) 449-517. */
    RP_ACN_TFA_COIL_BOLTZMANN, 
    //! Reversed phase, ACN, formic acid, ROD_BOLTZMANN model.
    /*! A ChemicalBasis calibrated for reversed phase, ACN as a second solvent,
        0.1% FA and ROD_BOLTZMANN type of BioLCCC model. The data was obtained
        in the joint research of Harvard University and Institute for Energy 
        Problems for Chemical Physics, Russian Academy of Science. */
    RP_ACN_FA_ROD_BOLTZMANN 
};

//! An instance of ChemicalBasis contains a set of BioLCCC constants.
/*!
    An instance of ChemicalBasis manages all the physicochemical constants,
    which are used in the calculations. Currently it contains:
        - The list of amino acids and peptide terminal groups.
        - The terminal groups which are set by default (cannon be changed).
        - The Energy of binding between a solvent and the surface of a solid 
          phase.
        - The type of BioLCCC model being used in calculations. 
        - Peptide geometry: the length of amino acid and the Kuhn length.
        - The width of the adsorbing layer.
       
    Note that the set of constants is highly interconnected. Usually the change
    in one constant, like the width of adsorbing layer or type of BioLCCC model,
    would only deteriorate the quality of RT prediction.
 */
class ChemicalBasis
{
public:
    //! Constructs an empty ChemicalBasis instance.
    ChemicalBasis();

    //! Constructs a ChemicalBasis instance with a predefined set of constants.
    ChemicalBasis(PredefinedChemicalBasis predefinedChemicalBasisId);

    //! Returns the map of all chemical groups.
    /*!
        A chemical group can be retrieved from the map by its label.
     */
    const std::map<std::string, ChemicalGroup> & chemicalGroups() const;

    //! Returns the default N-terminal group.
    const ChemicalGroup & defaultNTerminus() const;

    //! Returns the default C-terminal group.
    const ChemicalGroup & defaultCTerminus() const;

    //! Adds \a newChemicalGroup to the ChemicalBasis.
    /*!
        If an instance of ChemicalBasis already contains a chemical group with
        the same label than it is overwritten.
     */
    void addChemicalGroup(ChemicalGroup newChemicalGroup);

    //! Removes a chemical group with the given \a label;
    /*!
        Throws ChemicalBasisException if a chemical group with the given label
        is not found.
    */
    void removeChemicalGroup(std::string label);

    //! Removes all chemical groups in a ChemicalBasis.
    void clearChemicalGroups();

    //! Sets \a newBindEnergy as the binding energy of chemical group \a label.
    /*!
        Throws ChemicalBasisException if the chemical group is not found.
        \param label The label of the chemical group to be modified.
        \param newBindEnergy The new value of the bind energy.
    */
    void setChemicalGroupBindEnergy(std::string label, double newBindEnergy);

    //! Returns the bind energy of the second solvent. 
    /*! 
        Note that the bind energy of water is zero and the unit is kT.
    */
    double secondSolventBindEnergy() const;

    //! Sets \a newEnergy as the bind energy of the second solvent. 
    /*!
        Note that the bind energy of water is zero and the unit is kT.
    */
    void setSecondSolventBindEnergy(double newEnergy);

    //! Sets the type of BioLCCC model (e.g. CoilBoltzmann, CoilSnydel).
    void setModel(ModelType newModel);

    //! Returns the type of BioLCCC model (e.g. CoilBoltzmann, CoilSnydel).
    const ModelType model() const;

    //! Returns the length of a single amino acid residue in angstroms.
    /*!
        Due to the complex geometry of peptide molecule, this length is defined
        only approximately. The definition is the average length of an amino
        acid residue along backbone. In other terms, it is the length of 
        a backbone divided by the number of amino acid residues.
     */
    double segmentLength() const;

    //! Sets the length between two peptide bonds in angstroms.
    /*!
        Due to the complex geometry of peptide molecule, this length is defined
        only approximately. The definition is the average length of an amino
        acid residue along backbone. In other terms, it is the length of 
        a backbone divided by the number of amino acid residues.
     */
    void setSegmentLength(double newSegmentLength);

    //! Returns the Kuhn length of a molecule in amino acid residues.
    /*!
        A polymer molecule can be modelled as a chain of equal-sized rigid rods,
        freely joined with each other. In this case, the rods would be called 
        Kuhn segments and the length of a segment would be the Kuhn length.
        In BioLCCC model the Kuhn length is measured as a number of amino acid 
        residues comprising one rigid segment.

        This value is used only for COIL_* types of model.
    */
    int kuhnLength() const;

    //! Sets the Kuhn length of a molecule in amino acid residues.
    /*!
        A polymer molecule can be modelled as a chain of equal-sized rigid rods,
        freely joined with each other. In this case, the rods would be called 
        Kuhn segments and the length of a segment would be the Kuhn length.
        In BioLCCC model the Kuhn length is measured as a number of amino acid 
        residues comprising one rigid segment.

        This value is used only for COIL_* types of model.
    */
    void setKuhnLength(int newKuhnLength);

    //! Returns the width of a solid phase adsorption layer.
    /*!
        The width of a solid phase adsorption layer can be defined as a
        characteristic distance of interaction between an amino acid residue and
        the surface of a solid phase. 

        This value is used only for ROD_* types of model.
    */
    double adsorbtionLayerWidth() const;

    //! Sets the width of a solid phase adsorption layer.
    /*!
        The width of a solid phase adsorption layer can be defined as a
        characteristic distance of interaction between an amino acid residue and
        the surface of a solid phase. 

        This value is used only for ROD_* types of model.
    */
    void setAdsorbtionLayerWidth(double newAdsorbtionLayerWidth);

    //! Sets one of predefined chemical basis.
    ChemicalBasis setPredefinedChemicalBasis(
        PredefinedChemicalBasis predefinedChemicalBasisId);

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
