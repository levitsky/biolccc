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

The Python examples are specific to Python 2.x, since our project doesn't
support Python 3.x.

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
    BioLCCC model. The initial data were taken from Guo et al, Journal of 
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

Peptide sequence notation
*************************

In libBioLCCC we use the extended peptide notation. It is based on the
`one-letter IUPAC notation <http://www.chem.qmul.ac.uk/iupac/AminoAcid/>`_, 
but borrows only letters for the standard 20 aminoacid (i.e. no B, Z, X). 
We extended it in the following way:

- Modified amino acids are denoted as **xyzX**, i.e. their labels start with an 
  arbitrary number of lower-case letters and terminate with a single
  upper-case letter. The upper-case letter shows the base amino acid, while the
  lower-case letters describe the type of modification. The examples are:

    - **oxM** for oxidated methionine
    - **pS** for phosphorylated serine
    - **pT** for phosphorylated threonine
    - **camC** for carboxyamidomethylated cysteine

- The non-standard peptide terminal groups are denoted as **XxXx-** and
  **-XxXx**
  for N-terminal and C-terminal groups correspondingly. The label could contain
  an arbitrary number of mixed lower-case and upper-case letters and numbers, 
  but it should not be
  a valid peptide sequence. If a terminal group is not specified, it is
  assumed to be the standard one (i.e. an N-terminal hydrogen atom or C-terminal
  acidic group). The examples:
  
    - **Ac-** for N-Terminal acetylation
    - **H-** for N-Terminal hydrogen
    - **-NH2** for C-Terminal amidation
    - **-OH** for C-Terminal carboxyl group

- If a sequence contains two dots, then only the substring between them is
  parsed. This notation is used in several MS/MS search engines to show the
  adjacent amino acid residues for a peptide cleaved out of a protein. The
  examples are:

    -  K.APGFGDNR.K
    -  K.VGEVIVTK.D

Calculating retention time
**************************

calculateRT is the first libBioLCCC function you may need.
It requires three arguments: a peptide sequence,
a chemical basis, and and a description of chromatographic conditions. Supplied 
with these data, it
calculates the retention time of the peptide.

.. list-table:: Calculating the retention time of a peptide
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. code-block:: cpp

          #include <iostream>
          #include <string>
          #include <biolccc.h>

          int main() {
              std::string peptide("Ac-PEPTIDE-NH2");
              double RT = BioLCCC::calculateRT(peptide,
                  BioLCCC::rpAcnFaRodBoltzmann,
                  BioLCCC::standardChromoConditions);
              std::cout << "The retention time of " 
                        << peptide << " is " << RT << std::endl;
              return 0;
          }

     - 

       .. code-block:: python

          import pyBioLCCC

          peptide = 'Ac-PEPTIDE-NH2'
          RT = pyBioLCCC.calculateRT(peptide,
                   pyBioLCCC.rpAcnFaRodBoltzmann,
                   pyBioLCCC.standardChromoConditions)
          print 'The retention time of', peptide, 'is', RT

Please, consult with the libBioLCCC documentation for the details of calculateRT
function.

Specifying chromatographic conditions
*************************************

The next thing you may need to learn is how to specify the chromatographic
conditions. In order to do that, create a new instance of ChromoConditions and
replace the default parameters with your own.

.. list-table:: Specifying chromatographic conditions
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. code-block:: cpp

          #include <iostream>
          #include <biolccc.h>

          int main() {
              myChromoConditions = BioLCCC::ChromoConditions()

              // The column length in mm.
              myChromoConditions.setColumnLength(100.0);

              // The internal column diameter in mm.
              myChromoConditions.setColumnDiameter(0.1);

              // The average pore size in A.
              myChromoConditions.setColumnPoreSize(300.0);

              // The concentration of the eluting solvent (ACN for the reversed
              // phase) in component A in %.
              myChromoConditions.setSecondSolventConcentrationA(5.0);

              // The concentration of the eluting solvent (ACN for the reversed
              // phase) in component B in %.
              myChromoConditions.setSecondSolventConcentrationB(80.0);

              // The shape of the gradient. The example is a linear gradient
              // from 0% to 90% of component B over 60 minutes.
              myChromoConditions.setGradient(
                  BioLCCC::Gradient(0.0, 90.0, 60.0));
              
              // The flow rate in ml/min. 
              myChromoConditions.setFlowRate(0.0005);

              std::string peptide("Ac-PEPTIDE-NH2");
              double RT = BioLCCC::calculateRT(peptide,
                  BioLCCC::rpAcnFaRodBoltzmann,
                  myChromoConditions);
              std::cout << "The retention time of " 
                        << peptide << " is " << RT << std::endl;
              return 0;
          }

     - 

       .. code-block:: python

          import pyBioLCCC

          myChromoConditions = pyBioLCCC.ChromoConditions()

          # The column length in mm.
          myChromoConditions.setColumnLength(100.0)

          # The internal column diameter in mm.
          myChromoConditions.setColumnDiameter(0.1)

          # The average pore size in A.
          myChromoConditions.setColumnPoreSize(300.0)

          # The concentration of the eluting solvent (ACN for the reversed
          # phase) in component A in %.
          myChromoConditions.setSecondSolventConcentrationA(5.0)

          # The concentration of the eluting solvent (ACN for the reversed
          # phase) in component B in %.
          myChromoConditions.setSecondSolventConcentrationB(80.0)

          # The shape of the gradient. The example is a linear gradient
          # from 0% to 90% of component B over 60 minutes.
          myChromoConditions.setGradient(pyBioLCCC.Gradient(0.0, 90.0, 60.0))

          # The flow rate in ml/min. 
          myChromoConditions.setFlowRate(0.0005)

          peptide = 'Ac-PEPTIDE-NH2'
          RT = pyBioLCCC.calculateRT(peptide,
                   pyBioLCCC.rpAcnFaRodBoltzmann,
                   myChromoConditions)
          print 'The retention time of', peptide, 'is', RT

pyBioLCCC adds another way to interact with ChromoConditions. You can use its
instances as Python dictionaries:

.. list-table:: Dict-like syntax of ChromoConditions
   :widths: 40
   :header-rows: 1

   * - Python
   * - 

       .. code-block:: python

          import pyBioLCCC

          myChromoConditions = pyBioLCCC.ChromoConditions()
          print myChromoConditions.keys()

          myChromoConditions['columnLength'] = 100.0
          myChromoConditions['columnDiameter'] = 0.1
          myChromoConditions['columnPoreSize'] = 300.0
          myChromoConditions['secondSolventConcentrationA'] = 5.0
          myChromoConditions['secondSolventConcentrationB'] = 80.0
          myChromoConditions['gradient'] = pyBioLCCC.Gradient(0.0, 90.0, 60.0)
          myChromoConditions['flowRate'] = 0.0005

          peptide = 'Ac-PEPTIDE-NH2'
          RT = pyBioLCCC.calculateRT(peptide,
                   pyBioLCCC.rpAcnFaRodBoltzmann,
                   myChromoConditions)
          print 'The retention time of', peptide, 'is', RT

Besides being more convenient and compact, this syntax allows ChromoConditions 
to be pickled. 

If you want to see the full list of parameters stored in a ChromoConditions
instance, please, take a look at the class description in the libBioLCCC
documentation.

Calculating mass
****************

libBioLCCC contains functions to calculate the monoisotopic and average masses
of a peptide. Besides the sequence of a peptide, you need to supply a
ChemicalBasis instance which contains the masses of amino acids. 

.. list-table:: Calculating mass of a peptide
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. code-block:: cpp

          #include <iostream>
          #include <biolccc.h>

          int main() {

              std::string peptide("Ac-PEPTIDE-NH2");

              double averageMass = BioLCCC::calculateAverageMass(
                  peptide, BioLCCC::rpAcnFaRodBoltzmann);
              double monoisotopicMass = BioLCCC::calculateMonoisotopicMass(
                  peptide, BioLCCC::rpAcnFaRodBoltzmann);

              std::cout << "Average mass of " << peptide << " is " 
                        << averageMass << " Da" << std::endl;
              std::cout << "Monoisotopic mass of " << peptide << " is " 
                        << monoisotopicMass << " Da" << std::endl;

              return 0;
          }

     - 

       .. code-block:: python

          import pyBioLCCC

          peptide = 'Ac-PEPTIDE-NH2'

          averageMass = pyBioLCCC.calculateAverageMass(
              peptide, pyBioLCCC.rpAcnFaRodBoltzmann)
          monoisotopicMass = pyBioLCCC.calculateMonoisotopicMass(
              peptide, pyBioLCCC.rpAcnFaRodBoltzmann)

          print 'The average mass of', peptide, 'is', averageMass, 'Da'
          print 'The monoisotopic mass of',peptide, 'is', monoisotopicMass, 'Da'


..
    .. list-table:: example of a code snippet
       :widths: 40 40
       :header-rows: 1

       * - C++
         - Python
       * - 

           .. code-block:: cpp

              #include <iostream>

              int main() {
                  return 0;
              }

         - 

           .. code-block:: python

              print a+b

