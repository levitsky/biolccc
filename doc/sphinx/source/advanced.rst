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

       .. literalinclude:: ../../../src/examples/kd_calculation.cpp
          :language: cpp

     - 

       .. literalinclude:: ../../../src/examples/kd_calculation.py
          :language: python


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

       .. literalinclude:: ../../../src/examples/gradient.cpp
          :language: cpp

     - 

       .. literalinclude:: ../../../src/examples/gradient.py
          :language: python

Changing chemical bases
***********************

If you need to change the predefined values of physicochemical constants, you
may edit an instance of ChemicalBasis. The 
`ChemicalBasis class documentation <./API/classBioLCCC_1_1ChemicalBasis.html>`_
contains further information on all the parameters it contains.

Please, keep in mind a few rules on editing a ChemicalBasis instance:

- **You cannot change the label of a ChemicalGroup.**

  This is happening, because the label of a ChemicalGroup is also used in a
  ChemicalGroups map of a ChemicalBasis. If you still want for some reason 
  change modify it, you need to remove this group and add a new one with the
  same chemical properties but another label.

- **Do not modify the predefined chemical bases.**
 
  You cannot broke this rule in C++, since the predefined bases made constant
  there. But in Python there are no constants, and you can accidentally 
  modify variables pyBioLCCC.rpAcnFaRod and 
  pyBioLCCC.rpAcnTfaCoil. You should avoid this because you can easily
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
  BioLCCC model, you need to conduct the whole calibration procedure which
  includes the LC-experiment with the calibration mixture and the further data
  processing.

Here is an example of code modifying a ChemicalBasis instance:

.. list-table:: Modifying a ChemicalBasis instance
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: ../../../src/examples/advanced_chembasis.cpp
          :language: cpp

     - 

       .. literalinclude:: ../../../src/examples/advanced_chembasis.py
          :language: python


Once again, the libBioLCCC documentation contains the full list of available
parameters of ChemicalBasis.

Parsing peptide sequence
************************

Sequence parsing is a process in which a text sequence is translated into a list
of chemical groups. The resulting list begins with the N-terminal group,
continues with the amino acids and ends with the C-Terminal group.

.. list-table:: Parsing peptide sequence
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: ../../../src/examples/sequence_parsing.cpp
          :language: cpp

     - 

       .. literalinclude:: ../../../src/examples/sequence_parsing.py
          :language: python


Changing the calculation precision
**********************************

The main equation of liquid chromatography involves the intergration over the
pumped volume of binary solvent:

.. math::

   \int_{0}^{V_R - V_0}{\frac{dV}{V_P \, K_D(V)}} = 1

where *V* is the volume of binary solvent pumped through the column, 
*V*\ :sub:`R` is the
retention volume of a substance, *V*\ :sub:`P` is the volume of pores and 
*V*\ :sub:`0` is the dead volume of the chromatographic system.

libBioLCCC computes this integral as a summation over values of V. The step of
this summation is dV. You can change this value using an instance of
ChromoConditions. By default, dV equals zero, which means that its value is
derived from the flow rate. Currently, if dV == 0 than 
:math:`dV = flow\; rate/20`

Advanced customization
**********************

Segmentation mechanism
======================

BioLCCC uses two type of units to divide a polymer molecule. The first is a
conventional monomer, i.e. a building block of a molecule. In case of proteins
it is an amino acid residue. The terminal groups are not considered as monomers,
rather they are modifiers attached to monomers.

But when we want to describe the conformations of a molecule, a monomer is not
always a good unit. The standard model of a long polymer molecule is a chain of
free jointed rigid rods, or Kuhn segments. The length of a Kuhn segment does not
necessarily equal to the length of a monomer, and it even may not be its 
multiple.
That is why we use in BioLCCC an additional scale, the Kuhn length of a polymer.
The Kuhn length is a minimal distance between two chemical bonds in a polymer
backbone, whose orientations are *almost* independent of each other.

The calculation procedure is the following. At first, we define the sequence of
monomers in a polymer chain and calculate their effective adsorption energies. 
Then we divide the chain into Kuhn segments and assign each monomers to the
corresponding segment. If a boundary between Kuhn segments crosses the monomer
then the monomer itself is divided into two parts, and each is assigned to the
corresponding segment. The effective adsorption of a Kuhn segment is a sum of
effective energies of monomers it contains. If a segment contains only a part of
monomer then its energy is taken proportional to its length. 

In the case of COIL model, the centers of this segments are modelled as
adsorbing beads which are connected by freely jointed rods. 
For the ROD model, the centers of segments become the beads threaded regularly 
on a single rigid rod. The distance between the beads in both cases equals to
the Kuhn length.

Adsorption factors in the COIL model
====================================

The standard COIL BioLCCC model assumes that adsorption occurs only in
a single layer, which is closest to the wall. In this case
This assumption can be generalized to the case when several near-wall layers
adsorb the segments of a polymer chain. In terms of translational matrices it
means, that the second and further rows would contain the exponential terms.
Because the energy of binding to these layers may differ, we introduce
the vector of layer-specific values of adsorption strength. It is contained in 
adsorptionLayerFactors() function of a ChemicalBasis instance. 
The first element of the vector corresponds to the layer closest to the wall, 
second to the next and so on. The vector may contain an arbitrary number of
elements, but it must be less than a half of the number of rows in the
transitional matrix. This number is calculated as (pore size / kuhn length).

Note that the relative adsorption strength of the column is applied to these
values as usual.

Other custom options
====================

- You may enable the linear approximation of binary solvent effective energy
  using the function setSnyderApproximation(true).
- If you are using the solvents other than water/ACN, you may set the
  corresponding densities and molar masses using the functions of ChemicalBasis.
  These values are used in the equation for the effective binding energy of
  a binary solvent.

