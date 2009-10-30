#ifndef CHROMOCONDITION_H
#define CHROMOCONDITION_H

#include <string>
#include <utility>
#include <vector>

#include "gradient.h"

//! This class encapsulates all parameters of chromatographic equipment.

class ChromoConditions {  
    public:
        
        /*!
            Constructs a custom ChromoConditions object with the standard 
            Dionex gradient conditions.
        */
        ChromoConditions(double iColumnLength = 150.0,
                         double iColumnDiameter = 0.075,
                         double iColumnPoreSize = 100.0,
                         Gradient iGradient = Gradient(0.0, 50.0, 60.0),
                         double iSecondSolventConcentrationA = 2.0,
                         double iSecondSolventConcentrationB = 80.0,
                         double iDelayTime = 0.0,
                         std::string iSecondSolvent = std::string("ACN"),
                         double iFlowRate = 0.0003,
                         double iDV = 0.0,
                         double iCalibrationParameter = 1.0,
                         double iColumnVpToVtot = 0.5,
                         double iColumnPorosity = 0.9,
                         double iTemperature = 293.0);
        
        /*!
            Returns the length of the column in mm.
        */
        double columnLength() const;
        
        /*!
            Sets the length of the column in mm.
        */
        void setColumnLength(double newColumnLength);
        
        /*!
            Returns the internal diameter of the column in mm.
        */
        double columnDiameter() const;
        
        /*!
            Sets the internal diameter of the column in mm.
        */
        void setColumnDiameter(double newColumnDiameter);
        
        /*!
            Returns the size of the pores in angstroms.
        */
        double columnPoreSize() const;
        
        /*!
            Sets the size of the pores in angstroms.
        */
        void setColumnPoreSize(double newColumnPoreSize);
        
        /*!
            Returns the volume of pores divided by the total column volume 
            (Pi * r^2 * l).
        */
        double columnVpToVtot() const;
        
        /*!
            Sets the volume of the pores divided by the total column volume 
            (Pi * r^2 * l).
        */
        void setColumnVpToVtot(double newColumnVpToVtot);
        
        /*!
            Returns ( volume of pores + volume of liquid phase ) / 
            total column volume ( PI * r^2 * l)
        */
        double columnPorosity() const;
        
        /*!
            Sets ( volume of pores + volume of liquid phase ) / 
            total column volume ( PI * r^2 * l)
        */
        void setColumnPorosity(double newColumnPorosity);
        
        /*!
            Returns the temperature in degrees of Kelvin.
        */
        double temperature() const;
        
        /*!
            Sets the temperature in degrees of Kelvin.
        */
        void setTemperature(double newTemperature);
        
        /*!
            Returns the calibration parameter.
        */
        double calibrationParameter() const;
        
        /*!
            Sets the calibration parameter.
        */
        void setCalibrationParameter(double newCalibrationParameter);
        
        /*!
            Returns the flow rate in ml/min.
        */
        double flowRate() const;
        
        /*!
            Sets the flow rate in ml/min.
        */
        void setFlowRate(double newFlowRate);
        
        /*!
            Returns the volume of the pump mixer, ml.
        */
        double dV() const;
        
        /*!
            Sets the volume of the pump mixer, ml.
        */
        void setDV(double newDV);
        
        /*!
            Returns the delay time;
        */
        double delayTime() const;
        
        /*!
            Sets the delay time;
        */
        void setDelayTime(double newDelayTime);
        
        /*!
            Returns the name of the second solvent.
        */
        std::string secondSolvent() const;
        
        /*!
            Sets the name of the second solvent.
        */
        //std::string setSecondSolvent(std::string);
        
        /*!
            Returns the concentration of the second solvent in the component A.
        */
        double secondSolventConcentrationA() const;

        /*!
            Sets the concentration of the second solvent in the component A.
        */
        void setSecondSolventConcentrationA(
            double newSecondSolventConcentrationA);

        /*!
            Returns the concentration of the second solvent in the component B.
        */
        double secondSolventConcentrationB() const;

        /*!
            Sets the concentration of the second solvent in the component B.
        */
        void setSecondSolventConcentrationB(
            double newSecondSolventConcentrationB);
        
        /*!
            Returns the elution gradient.
        */
        Gradient gradient() const;

        /*!
            Sets the elution gradient.
        */
        void setGradient(Gradient newGradient);
        
    private:
        double mColumnLength;
        double mColumnDiameter;
        double mColumnPoreSize;
        double mColumnVpToVtot;
        double mColumnPorosity;
        double mTemperature;
        double mCalibrationParameter;
        double mFlowRate;
        double mDV;
        double mDelayTime;
        Gradient mGradient;
        std::string mSecondSolvent;
        double mSecondSolventConcentrationA;
        double mSecondSolventConcentrationB;
};

#endif
