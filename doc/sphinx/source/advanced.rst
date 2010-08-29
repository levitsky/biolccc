=========================
Advanced libBioLCCC usage
=========================

Calculating the coefficient of distribution
*******************************************

libBioLCCC allows to calculate the coefficient of distribution of a peptide.
It is a fundamental measure of adsorption which is expressed as:

.. math::

   K_d = \frac{C_{peptide \: in\: pores}}
   {C_{peptide \: in \: intestitial \: volume}}

where *C* denotes concentration.

The function calculateKd does the calculations and requires the concentration of
the second solvent, a chemical basis, and the size of the pores. The optional
arguments are the relative strength of a column and temperature, but keep in
mind that they are experimental.

.. list-table:: Calculating the distribution coefficient
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

              double kd = BioLCCC::calculateKd(
                  peptide, // the peptide sequence 
                  15.0,    // the concentration of the second solvent, %
                  BioLCCC::rpAcnFaRodBoltzmann, // the chemical basis
                  100.0);  // the size of the pores, angstroms

              std::cout << "The coefficient of distribution of " << peptide << 
                        << " is " << kd << std::endl;
              return 0;
          }

     - 

       .. code-block:: python

          import pyBioLCCC

          peptide = 'Ac-PEPTIDE-NH2'

          kd = pyBioLCCC.calculateKd(
              peptide, # the peptide sequence 
              15.0,    # the concentration of the second solvent, %
              pyBioLCCC.rpAcnFaRodBoltzmann, # the chemical basis
              100.0)   # the size of the pores, angstroms

          print 'The coefficient of distribution of', peptide, 'is', kd

Non-linear gradients
********************

The shape of a gradient describes how the concentration of component B changes
over time. It is implemented as a list of points, each point telling what the
concentration of component B should be at given time. The solvent composition
between these points is calculated using linear interpolation.

.. list-table:: 
   :widths: 40 40
   :header-rows: 0

   * - For example, we want to describe a gradient consisting of four parts:

       - isocratic elution at 5% of component B during the first 20 minutes
       - a linear part raising from 5% to 45% during 40 minutes
       - a sharp rise to 100% of component B during 5 minutes
       - the final flushing step at 100% of component B over 20 minutes.

     -
      .. image:: ./gradient.png
         :scale: 50 %

The following code does that:

.. list-table:: Using non-linear gradients
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
              myChromoConditions = BioLCCC::ChromoConditions();
              myGradient = BioLCCC::Gradient();
              myGradient.addPoint(0.0, 5.0);
              myGradient.addPoint(20.0, 5.0);
              myGradient.addPoint(60.0, 45.0);
              myGradient.addPoint(65.0, 100.0);
              myGradient.addPoint(85.0, 100.0);
              myChromoConditions.setGradient(myGradient);

              double RT = BioLCCC::calculateRT(peptide,
                  BioLCCC::rpAcnFaRodBoltzmann,
                  myChromoConditions);
              std::cout << "The retention time of " 
                        << peptide << " in the custom gradient is " 
                        << RT << std::endl;
              return 0;
          }

     - 

       .. code-block:: python

          import pyBioLCCC

          peptide = 'Ac-PEPTIDE-NH2'
          myChromoConditions = pyBioLCCC.ChromoConditions()
          myGradient = pyBioLCCC.Gradient()
          myGradient.addPoint(0.0, 5.0)
          myGradient.addPoint(20.0, 5.0)
          myGradient.addPoint(60.0, 45.0)
          myGradient.addPoint(65.0, 100.0)
          myGradient.addPoint(85.0, 100.0)
          myChromoConditions.setGradient(myGradient)

          RT = pyBioLCCC.calculateRT(peptide,
                   pyBioLCCC.rpAcnFaRodBoltzmann,
                   myChromoConditions)
          print 'The retention time of', peptide, 'in the custom gradient is',RT

Changing chemical bases
***********************

If you need to change the predefined values of physicochemical constants, you
may edit an instance of ChemicalBasis. The ChemicalBasis class documentation
contains further information on all the parameters it contains.

However, while editing ChemicalBasis please keep in mind the few several rules:

- **Do not change the label of a chemical group.**

  The function .setLabel() of
  ChemicalGroup is made only for technical reasons and should not be invoked by
  a user. The same goes for Python dict-like syntax which may somehow change the
  label of a ChemicalGroup (for example, 
  pyBioLCCC.rpAcnFaRodBoltzmann['A']['label'] = 'B'). If you need to change the
  label of a chemical group, please create a new ChemicalGroup instance with the
  parameters you need.

- **Do not modify the predefined chemical bases.**
 
  You cannot broke this rule in C++, since the predefined bases made constant
  there. However, in Python there are no constants and you can accidentally 
  modify variables pyBioLCCC.rpAcnFaRodBoltzmann and 
  pyBioLCCC.rpAcnTfaCoilBoltzmann. You should avoid this because you can easily
  forget about it later and use these bases as if they were intact.

  If you need to derive a new basis from a predefined one, use an alternative
  ChemicalBasis constructor. This constructor requires a name of predefined
  chemical basis and fills a newly created instance with the corresponding data.

  The names of predefined chemical bases contain in PredefinedChemicalBasis
  type. For the further information, please consult libBioLCCC documentation.

- **Arbitrary changes in a chemical basis are likely to worsen the accuracy of 
  prediction.**

  The constants stored in a chemical basis were found using a combination
  of specially developed LC experiments and calibration algorithms. These data
  correspond to the local maximum of predicting ability for a given combination
  of solvents and stationary phase. According to our experience, it is unlikely
  that a change in a single constant will rise the accuracy of prediction. If
  you need to adopt BioLCCC to a custom retention chemistry or another type of
  BioLCCC model, you need to measure the 
  retention times of the calibration mixture and recalculate the ChemicalBasis.

Here is an example of code modifying a ChemicalBasis instance:

.. list-table:: Modifying a ChemicalBasis
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. code-block:: cpp

          #include <iostream>
          #include <biolccc.h>

          int main() {
              // Deriving a new ChemicalBasis instance from a predefined one.
              BioLCCC::ChemicalBasis
                  myChemicalBasis(BioLCCC::RP_ACN_FA_ROD_BOLTZMANN);

              // Changing the bind energy of a chemical group.
              myChemicalBasis.chemicalGroups()["E"].setBindEnergy(0.0);
              myChemicalBasis.chemicalGroups()["-NH2"].setBindEnergy(0.0);

              std::cout << "The bind energy of E is "
                  << myChemicalBasis.chemicalGroups()["E"].bindEnergy()
                  << std::endl;
              std::cout << "The bind energy of -NH2 is "
                  << myChemicalBasis.chemicalGroups()["-NH2"].bindEnergy()
                  << std::endl;

              // Adding a new chemical group. The energy is not valid.
              myChemicalBasis.addChemicalGroup(
                  BioLCCC::ChemicalGroup(
                      "Hydroxyproline",      // full name
                      "hoP",                 // label
                      0.40,                  // bind energy
                      97.1167+15.9994,       // average mass
                      97.05276+15.9994915)); // monoisotopic mass

              // Setting a new type of model. Without a massive recalibration
              // it will ruin the accuracy of prediction.
              myChemicalBasis.setModel(BioLCCC::COIL_BOLTZMANN);

              std::string peptide("Ac-PEhoPTIDE-NH2");
              double RT = BioLCCC::calculateRT(peptide,
                  myChemicalBasis,
                  BioLCCC::standardChromoConditions);

              double monoisotopicMass = BioLCCC::calculateMonoisotopicMass(
                  peptide, myChemicalBasis);

              std::cout << "The retention time of " 
                        << peptide << " is " << RT << std::endl;
              std::cout << "The monoisotopic mass of " << peptide << " is " 
                        << monoisotopicMass << " Da" << std::endl;

              return 0;
          }

     - 

       .. code-block:: python

          import pyBioLCCC

          # Deriving a new ChemicalBasis instance from a predefined one.
          myChemicalBasis pyBioLCCC.ChemicalBasis(
              pyBioLCCC.RP_ACN_FA_ROD_BOLTZMANN)

          # Changing the bind energy of a chemical group.
          myChemicalBasis.chemicalGroups()['E'].setBindEnergy(0.0)
          myChemicalBasis.chemicalGroups()['-NH2'].setBindEnergy(0.0)

          print "The bind energy of E is", \
              myChemicalBasis.chemicalGroups()['E'].bindEnergy()
          print "The bind energy of -NH2 is", \
              myChemicalBasis.chemicalGroups()['-NH2'].bindEnergy()

          # Adding a new chemical group. The energy is not valid.
          myChemicalBasis.addChemicalGroup(
              pyBioLCCC.ChemicalGroup(
                  'Hydroxyproline',      # full name
                  'hoP',                 # label
                  0.40,                  # bind energy
                  97.1167+15.9994,       # average mass
                  97.05276+15.9994915))  # monoisotopic mass

          # Setting a new type of model. Without a massive recalibration
          # it will ruin the accuracy of prediction.
          myChemicalBasis.setModel(pyBioLCCC.COIL_BOLTZMANN);

          peptide = "Ac-PEhoPTIDE-NH2"
          RT = pyBioLCCC.calculateRT(peptide,
              myChemicalBasis,
              pyBioLCCC.standardChromoConditions)

          monoisotopicMass = pyBioLCCC.calculateMonoisotopicMass(
              peptide, myChemicalBasis)

          print 'The retention time of', peptide, 'is', RT
          print 'The monoisotopic mass of', peptide, 'is', monoisotopicMass,'Da'

Parsing peptide sequence
************************

Changing the calculation precision
**********************************


