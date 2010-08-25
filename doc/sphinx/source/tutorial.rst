========
Tutorial
========

Before we begin
***************

The following help is written both for libBioLCCC and pyBioLCCC. The only
difference between two these packages lies in the syntax of commands. That is
why we supply code snippets both for C++ and Python. Here is an example:

.. list-table:: Example of a code snippet
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. code-block:: cpp

          #include <iostream>

          int main() {
              int a = 1;
              int b = 2;

              std::cout << a+b << std::endl;

              return 0;
          }

     - 

       .. code-block:: python

          a = 1
          b = 2

          print a+b

As you can see, the difference is noticeable, but it will be much smaller
for libBioLCCC examples. Note that we use Python 2.x, since Python 3.x is
still not widely adopted.

Basis conceptions
*****************

There are a few simple conceptions which are widely used in libBioLCCC. Most of
them are represented by a corresponding class. Here they are:

**BioLCCC model type** - a set of assumptions on the properties of protein
molecules. Strictly speaking, BioLCCC is a family of models, each based on the
different assumptions. This version of the BioLCCC model contains two types of
model:

    **ROD_BOLTZMANN** - in this type of model a peptide is represented as an
    absolutely rigid rod. Amino acids are modelled as regularly spaced beads
    threaded on this rod. This assumption works better for peptides rather
    than protein molecules. The adsorption is described by the Boltzmann
    equation.

    The equations for ROD_BOLZMANN are going to be published in the upcoming
    paper.
    
    **COIL_BOLTZMANN** - in this model a protein molecule is described as
    a flexible polymer. The conformations of this molecule in a pore can be
    modelled as a random walk in the field of adsorbing walls. This assumption
    should work better for long protein molecules. The adsorption itself is
    described by the Boltzmann equation.

    COIL_BOLTZMANN model was described in ''Liquid Chromatography at Critical 
    Conditions: Comprehensive Approach to Sequence-Dependent Retention Time 
    Prediction'', Alexander V. Gorshkov et al, Analytical Chemistry, 2006, 78
    (22), 7770-7777. `Link <http://dx.doi.org/10.1021/ac060913x>`_.

**Chemical group** - in libBioLCCC that is an amino acid residue OR a peptide
terminal group in a peptide chain. Examples are a histidine residue, 
phosphoserine residue and N-Terminal hydrogen that closes a peptide chain. The
properties of a chemical group are stored in the ChemicalGroup class. 

**Chemical basis** - a set of all physicochemical constants involved into the
BioLCCC equations. This set contains:

    - The list of all chemical groups, i.e. amino acids and terminal groups. 
      Any peptide can be represented as a series of these, that is why it is
      a *basis* similar to the mathematical basis. 
    - Which terminal groups are set by default.
    - The energy of binding between a solvent and the surface of a solid phase.
    - The type of BioLCCC model being used in calculations.
    - Peptide geometry: the length of an amino acid and the Kuhn length.
    - The range of an interaction between an amino acid and the surface of the
      solid phase (a.k.a. the width of the adsorbing layer).

The properties of a chemical basis are stored in the ChemicalBasis class.

A chemical basis is specific to a type of retention chemistry, solvents
and ion paring agent being used in the experiment. In addition, it must be used
only with the same type of model as the one used in the calibration of the
chemical basis.

**Predefined chemical basis** - a chemical basis, calculated (or, more
precisely, calibrated) for the specific retention chemistry and type of
BioLCCC model. The current version of libBioLCCC contains two predefined
chemical bases:

    **rpAcnFaRodBoltzmann** - a ChemicalBasis calibrated for a reversed phase,
    ACN as a second solvent, 0.1% FA in both solvents and ROD_BOLTZMANN type of
    BioLCCC model. The data was obtained in the joint research of Harvard 
    University and Institute for Energy Problems for Chemical Physics, 
    Russian Academy of Science.

    **rpAcnTfaCoilBoltzmann** - a chemical basis calibrated for a reversed
    phase,
    ACN as a second solvent, 0.1% TFA in both solvents and COIL_BOLTZMANN  
    BioLCCC model. The initial data was taken from Guo et al, Journal of 
    Chromatography, 359 (1986) 449-517.


**Chromatographic conditions** - a description of a chromatographic equipment 
and its settings. Contains:

    - The geometry of the column.
    - The properties of the adsorbent: average size of the pores, porosity
      (i.e. percentage of volume not filled with the solid phase),
      (volume of pores)/(total volume of column) ratio, relative adsorption
      strength.
    - Elution parameters: the shape of the gradient, the composition of
      components, flow rate, delay time.
    - The step of integration over volume.
    - Temperature of a column (EXPERIMENTAL).

The default values were set rather arbitrarily.


