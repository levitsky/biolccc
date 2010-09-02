#ifndef CHEMICALBASIS_H
#define CHEMICALBASIS_H

#include <map>
#include <vector>
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
        flexibility of a protein molecule. It use actively matrix calculations.
        */
    COIL,
    /*! The BioLCCC model with the assumption of absolute rigidity of 
        a protein molecule. It works better for short molecules. It use explicit
        expressions for Kd. */
    ROD
};

//! This enum describes the predefined sets of physicochemical constants.
/*!
    The BioLCCC library contains several predifined sets of physicochemical
    constants. Please note that usually changing only one parameter in a whole
    set of constants deteriorate the quality of prediction.
 */
enum PredefinedChemicalBasis
{
    //! Reversed phase, ACN, trifluoracetic acid, COIL model.
    /*! A ChemicalBasis calibrated for reversed phase, ACN as a second solvent,
        0.1% TFA and COIL type of BioLCCC model. The data was 
        obtained in Guo et al, Journal of Chromatography, 359 (1986) 449-517. */
    RP_ACN_TFA_COIL, 
    //! Reversed phase, ACN, formic acid, ROD model.
    /*! A ChemicalBasis calibrated for reversed phase, ACN as a second solvent,
        0.1% FA and ROD type of BioLCCC model. The data was obtained
        in the joint research of Harvard University and Institute for Energy 
        Problems for Chemical Physics, Russian Academy of Science. */
    RP_ACN_FA_ROD 
};

//! An instance of ChemicalBasis contains a set of BioLCCC constants.
/*!
    An instance of ChemicalBasis manages all the physicochemical constants,
    which are used in the calculations. Currently it contains:
        - the list of amino acids and peptide terminal groups;
        - the terminal groups which are set by default (cannon be changed);
        - the energy of binding between a solvent and the surface of a solid 
          phase;
        - the type of BioLCCC model being used in calculations and
          approximations
          used in
        - peptide geometry: the length of amino acid and the Kuhn length;
        - the width of the adsorbing layer.
       
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

    //! Returns the map of all chemical groups. Non-constant version.
    /*!
        A chemical group can be retrieved from the map by its label.
     */
    std::map<std::string, ChemicalGroup> & chemicalGroups();

    //! Returns the map of all chemical groups. Constant version.
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

    ////!Sets \a newBindEnergy as the binding energy of chemical group \a label.
    ///*!
    //    Throws ChemicalBasisException if the chemical group is not found.
    //    \param label The label of the chemical group to be modified.
    //    \param newBindEnergy The new value of the bind energy.
    //*/
    //void setChemicalGroupBindEnergy(std::string label, double newBindEnergy);

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

    //! Returns the length of a single monomer in a polymer chain in angstroms.
    /*!
        Due to the complex geometry of peptide molecule, this length is defined
        only approximately. The definition is the average length of an amino
        acid residue along backbone. In other terms, it is the length of 
        a backbone divided by the number of amino acid residues.
     */
    double monomerLength() const;

    //! Sets the length of a single monomer in a polymer chain in angstroms.
    /*!
        Due to the complex geometry of peptide molecule, this length is defined
        only approximately. The definition is the average length of an amino
        acid residue along backbone. In other terms, it is the length of 
        a backbone divided by the number of amino acid residues.
     */
    void setMonomerLength(double newMonomerLength);

    //! Returns the Kuhn length of a polymer molecule in angstroms.
    /*!
        A polymer molecule can be modelled as a chain of equal-sized rigid rods,
        freely joined with each other. In this case, the rods would be called 
        Kuhn segments, and the length of a segment would be the Kuhn length.
        The effective adsorption energy of a Kuhn segment equals to the total
        adsorption energy of all monomers that contains in this segment. If only
        a part of monomer contains in a segment than its energy is taken
        proportionally.

        There are no joints in the ROD model, the whole molecule is assumed to
        be shorter than a single Kuhn segment. However, in ROD model kuhnLength
        is still used to calculate the energy profile of a rod. The whole rod is
        divided into segments of kuhnLength and each segment transforms into an
        adsorbing bead. The effective energy of adsorption equals to the total
        effective energy of a segment, with the same expression as in the COIL
        model.
    */
    double kuhnLength() const;

    //! Sets the Kuhn length of a molecule in angstroms.
    /*!
        A polymer molecule can be modelled as a chain of equal-sized rigid rods,
        freely joined with each other. In this case, the rods would be called 
        Kuhn segments, and the length of a segment would be the Kuhn length.
        The effective adsorption energy of a Kuhn segment equals to the total
        adsorption energy of all monomers that contains in this segment. If only
        a part of monomer contains in a segment than its energy is taken
        proportionally.

        There are no joints in the ROD model, the whole molecule is assumed to
        be shorter than a single Kuhn segment. However, in ROD model kuhnLength
        is still used to calculate the energy profile of a rod. The whole rod is
        divided into segments of kuhnLength and each segment transforms into an
        adsorbing bead. The effective energy of adsorption equals to the total
        effective energy of a segment, with the same expression as in the COIL
        model.
    */
    void setKuhnLength(double newKuhnLength);

    //! Returns the width of a solid phase adsorption layer in ROD model.
    /*!
        The width of a solid phase adsorption layer can be defined as a
        characteristic distance of interaction between an amino acid residue and
        the surface of a solid phase. 

        This value is used only in the ROD model.
    */
    double adsorptionLayerWidth() const;

    //! Sets the width of a solid phase adsorption layer in ROD model.
    /*!
        The width of a solid phase adsorption layer can be defined as a
        characteristic distance of interaction between an amino acid residue and
        the surface of a solid phase. 

        This value is used only in the ROD model.
    */
    void setAdsorptionLayerWidth(double newAdsorptionLayerWidth);

    //! Returns the absorption factors of the near-wall layers in COIL model.
    /*! 
        The standard COIL BioLCCC model assumes that adsorption occurs only in
        one layer located close to the wall. However, this assumption can be
        generalized to the case when several near-wall layers adsorb segments 
        of a polymer chain. This vector contains the relative adsorbtion
        strengths of near-wall layers. This adsorbtion strength have the same
        meaning as the relative adsorbtion strength of a column and multiplyed
        by it. The first element of the vector corresponds to the
        layer closest to the wall, second to the next and so on.

        This value is used only in the COIL model.
     */
    const std::vector<double> & adsorptionLayerFactors() const;

    //! Sets the absorption factors of the near-wall layers in COIL model.
    /*! 
        The standard COIL BioLCCC model assumes that adsorption occurs only in
        one layer located close to the wall. However, this assumption can be
        generalized to the case when several near-wall layers adsorb segments 
        of a polymer chain. This vector contains the relative adsorption
        strengths of near-wall layers. This adsorption strength have the same
        meaning as the relative adsorption strength of a column and multiplyed
        by it. The first element of the vector corresponds to the
        layer closest to the wall, second to the next and so on.

        This value is used only in the COIL model.
     */
    void setAdsorptionLayerFactors(
        std::vector<double> newAdsorptionLayerFactors);

    //! Returns true if the energy of binary solvent is linearly fitted.
    /*!
        If the value of snyderApproximation is true then the energy of binary
        solvent is expressed by:
        E_{ab} = secondSolventBindEnergy * Nb
     */
    bool snyderApproximation() const;

    //! Enables the linear approximation of the energy of binary solvent.
    /*!
        If the value of snyderApproximation is true then the energy of binary
        solvent is expressed by:
        E_{ab} = secondSolventBindEnergy * Nb
     */
    void setSnyderApproximation(bool flag);

    //! Returns the density of the first solvent in kg/m^3.
    double firstSolventDensity() const;

    //! Sets the density of the first solvent in kg/m^3.
    void setFirstSolventDensity(double newFirstSolventDensity);

    //! Returns the density of the second solvent in kg/m^3.
    double secondSolventDensity() const;

    //! Sets the density of the second solvent in kg/m^3.
    void setSecondSolventDensity(double newSecondSolventDensity);

    //! Returns the molecular mass of the first solvent in g/mol.
    double firstSolventAverageMass() const;

    //! Sets the molecular mass of the first solvent in g/mol.
    void setFirstSolventAverageMass(double newFirstSolventAverageMass);

    //! Returns the molecular mass of the second solvent in g/mol.
    double secondSolventAverageMass() const;

    //! Sets the molecular mass of the second solvent in g/mol.
    void setSecondSolventAverageMass(double newSecondSolventAverageMass);

    //! Sets one of predefined chemical basis.
    ChemicalBasis setPredefinedChemicalBasis(
        PredefinedChemicalBasis predefinedChemicalBasisId);

private:
    std::map<std::string,ChemicalGroup> mChemicalGroups;
    double mSecondSolventBindEnergy;
    double mMonomerLength;
    double mKuhnLength;
    double mAdsorptionLayerWidth;
    std::vector<double> mAdsorptionLayerFactors;
    ModelType mModel;
    double mFirstSolventDensity;
    double mSecondSolventDensity;
    double mFirstSolventAverageMass;
    double mSecondSolventAverageMass;
    bool mSnyderApproximation;
};

}

#endif
