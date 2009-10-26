#ifndef PEPTIDEMETHODS_H
#define PEPTIDEMETHODS_H

#include "chemicalbasis.h"
#include "chromoconditions.h"

//! This class contains standard methods to calculate various peptide properties.
/*!
    <b> Known issues </b>
    - Due to fundamental ambiguity you couldn't parse properly the sequence for amidated histidine (H-NH2).
*/

class PeptideMethods {
    public:
        /*!
            Calculates a retention time for a given peptide sequence using given table of peptide chemicals and a chromatographic conditions.
        */
        static double calculateRTBioLCCC(const std::string &sequence,
                                        const ChemicalBasis &chemBasis,
                                        const ChromoConditions &conditions);
        
        /*!
            Calculates an average (molar) mass for a given peptide sequence using given table of peptide chemicals.
        */
        static double calculateAverageMass(const std::string &sequence,
                                        const ChemicalBasis &chemBasis);
                                
        /*!
            Calculates a monoisotopic mass for a given peptide sequence using given table of peptide chemicals.
        */
        static double calculateMonoisotopicMass(const std::string &sequence,
                                                const ChemicalBasis &chemBasis);
        /*!
            Calculates a distribution coefficient (Kd) for a given peptide sequence in a BioLCCC model using given table of peptide chemicals, a name of second solvent, its concentration, a size of adsorbent's pores, a calibration parameter and a temperature.
        */
        static double calculateKdBioLCCC (const std::string &sequence,
                                        const ChemicalBasis &chemBasis,
                                        const std::string &secondSolvent,
                                        const double secondSolventConcentration,
                                        const double columnPoreSize = 100.0,
                                        const double calibrationParameter = 1.0,
                                        const double temperature = 293.0);
        
        /*!
            Created as transient solution for fast calculation of RTBioLCCC and masses.
        */
        static bool calculatePeptideProperties(const std::string &sequence,
                                        const ChemicalBasis &chemBasis,
                                        const ChromoConditions &conditions,
                                        double *RTBioLCCC,
                                        double *averageMass,
                                        double *monoisotopicMass);


        
    private:
        static bool parseSequence(const std::string &source, 
                                const ChemicalBasis &chemBasis,
                                std::vector<Aminoacid> *parsedPeptideStructure,
                                Terminus *NTerminus,
                                Terminus *CTerminus,
                                std::vector<double> *peptideEnergyProfile = NULL);
        static double calculateRTBioLCCC (const std::vector<double> &peptideEnergyProfile,
                                const ChemicalBasis &chemBasis,
                                const ChromoConditions &conditions);
        static double calculateKdBioLCCC (const std::vector<double> &peptideEnergyProfile,
                                const ChemicalBasis &chemBasis,
                                const std::string secondSolvent,
                                const double secondSolventConcentration,
                                const double columnPoreSize = 100.0,
                                const double calibrationParameter = 1.0,
                                const double temperature = 293.0);
};

#endif
