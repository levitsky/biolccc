How to install pyteomics.biolccc?
---------------------------------

Linux (Debian/Ubuntu):

::

    sudo apt-get install python-setuptools python-dev
    sudo easy_install pip
    sudo pip install pyteomics.biolccc

Windows:

* Download pre-compiled binary packages from the
  `list <http://pypi.python.org/pypi/pyteomics.biolccc#downloads>`_.

  OR

* If you have Enthought Python Distribution / ActivePython, execute in the
  command line:

  ::

      easy_install pip
      pip install pyteomics.biolccc

What is BioLCCC?
----------------

BioLCCC (Liquid Chromatography of Biomacromolecules at Critical Conditions) is a
model describing the adsorption of protein molecules on porous media. Its
main application is retention time prediction in liquid chromatography, although
the list of potential applications can be easily extended. Contrary to the other
models of peptide/protein chromatography, BioLCCC starts from very basic
assumptions regarding flexibility of a polypeptide chain, the shape of a pore,
the type of interactions neglected, etc. Given these assumptions, the coefficient of
distribution (Kd) of a peptide between the solid and mobile phases can be
derived using the methods of statistical physics of macromolecules. Finally, the
retention time of a peptide is calculated from Kd using the basic equation of
gradient chromatography.

Owing to the physical basis of the BioLCCC model, it contains very few free
parameters. The retention properties of an amino acid are characterized by a
single number, which is essentially the energy of interaction between the amino
acid and the surface of solid phase in pure water+ion paring agent. Given this
small number of phenomenological parameters, the BioLCCC model can be easily
adapted for an arbitrary type of chromatography not limited by phase or solvent
types. Moreover, its extension to peptides with post-translational modifications
is straightforward as it was shown for the phosphorylated amino acids.

Several papers regarding BioLCCC model were published:

1. Liquid Chromatography at Critical Conditions:  Comprehensive Approach to
Sequence-Dependent Retention Time Prediction, Alexander V. Gorshkov, Irina A.
Tarasova, Victor V. Evreinov, Mikhail M. Savitski, Michael L. Nielsen, Roman A.
Zubarev, and Mikhail V. Gorshkov, Analytical Chemistry, 2006, 78 (22),
7770-7777. Link: http://dx.doi.org/10.1021/ac060913x.

2. Applicability of the critical chromatography concept to proteomics problems:
Dependence of retention time on the sequence of amino acids, Alexander V.
Gorshkov A., Victor V. Evreinov V., Irina A. Tarasova, Mikhail V. Gorshkov,
Polymer Science B, 2007, 49 (3-4), 93-107.
Link: http://dx.doi.org/10.1134/S1560090407030098.

3. Applicability of the critical chromatography concept to proteomics problems:
Experimental study of the dependence of peptide retention time on the sequence
of amino acids in the chain, Irina A. Tarasova, Alexander V. Gorshkov, Victor V.
Evreinov, Chris Adams, Roman A. Zubarev, and Mikhail V. Gorshkov, Polymer
Science A, 2008, 50 (3), 309.
Link: http://www.springerlink.com/content/gnh84v62w960747n/.

4. Retention time prediction using the model of liquid chromatography of
biomacromolecules at critical conditions in LC-MS phosphopeptide analysis,
Tatiana Yu. Perlova, Anton A. Goloborodko, Yelena Margolin, Marina L.
Pridatchenko, Irina A. Tarasova, Alexander V. Gorshkov, Eugene Moskovets,
Alexander R. Ivanov and Mikhail V. Gorshkov, Accepted to Proteomics.
Link: http://dx.doi.org/10.1002/pmic.200900837.

What is pyteomics.biolccc?
--------------------------

pyteomics.biolccc is an open source library, which implements the BioLCCC model
in the combination of Python and C++ programming languages.
It performs most BioLCCC-related tasks, such as:

* predicts the retention time of peptides and proteins in given
  chromatographic conditions;
* predicts the adsorption properties of protein molecules, namely coefficient of
  distribution between mobile and solid phase;
* manages elution conditions and physicochemical constants;
* calculates masses of peptides and proteins.

What is libBioLCCC?
-------------------

libBioLCCC is the C++ layer of pyteomics.biolccc. libBioLCCC can be used
separately from the Python wrappings and has a clean and well-documented API.

Where can I find more information?
----------------------------------

The project documentation is hosted at
http://theorchromo.ru/docs

The source code of pyteomics.biolccc and underlying libBioLCCC C++ library is
open and hosted at https://github.com/levitsky/biolccc.
