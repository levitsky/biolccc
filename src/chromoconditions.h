#ifndef CHROMOCONDITION_H
#define CHROMOCONDITION_H

#include <string>
#include <utility>
#include <vector>

//! This class encapsulates all parameters of chromatographic equipment.
typedef std::pair<double,double>     gradientPoint;
typedef std::vector< gradientPoint > gradientProfile;

class ChromoConditions {  
    public:
        
        /*!
            Constructs a custom ChromoConditions object with the standard Dionex
            gradient conditions setted by default.
        */
        ChromoConditions(double iColumnLength = 150.0,
                        double iColumnDiameter = 0.075,
                        double iColumnPoreSize = 100.0,
                        gradientProfile iGradient = gradientProfile(),
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
            Returns a length of a column in mm.
        */
        double columnLength() const;
        
        /*!
            Sets a length of a column in mm.
        */
        void setColumnLength(double newColumnLength);
        
        /*!
            Returns an internal diameter of a column in mm.
        */
        double columnDiameter() const;
        
        /*!
            Sets an internal diameter of a column in mm.
        */
        void setColumnDiameter(double newColumnDiameter);
        
        /*!
            Returns a size of pores in angstroms.
        */
        double columnPoreSize() const;
        
        /*!
            Sets a size of pores in angstroms.
        */
        void setColumnPoreSize(double newColumnPoreSize);
        
        /*!
            Returns a volume of pores divided by a total column volume 
            (Pi * r^2 * l).
        */
        double columnVpToVtot() const;
        
        /*!
            Sets a volume of pores divided by a total column volume 
            (Pi * r^2 * l).
        */
        void setColumnVpToVtot(double newColumnVpToVtot);
        
        /*!
            Returns a ( volume of pores + volume of liquid phase ) / 
            total column volume ( PI * r^2 * l)
        */
        double columnPorosity() const;
        
        /*!
            Sets a ( volume of pores + volume of liquid phase ) / 
            total column volume ( PI * r^2 * l)
        */
        void setColumnPorosity(double newColumnPorosity);
        
        /*!
            Returns a temperature in kelvin degrees.
        */
        double temperature() const;
        
        /*!
            Sets a temperature in kelvin degrees.
        */
        void setTemperature(double newTemperature);
        
        /*!
            Returns a calibration parameter.
        */
        double calibrationParameter() const;
        
        /*!
            Sets a calibration parameter.
        */
        void setCalibrationParameter(double newCalibrationParameter);
        
        /*!
            Returns a flow rate in ml/min.
        */
        double flowRate() const;
        
        /*!
            Sets a flow rate in ml/min.
        */
        void setFlowRate(double newFlowRate);
        
        /*!
            Returns a volume of a pump mixer, ml.
        */
        double dV() const;
        
        /*!
            Sets a volume of a pump mixer, ml.
        */
        void setDV(double newDV);
        
        /*!
            Returns a delay time;
        */
        double delayTime() const;
        
        /*!
            Sets a delay time;
        */
        void setDelayTime(double newDelayTime);
        
        /*!
            Returns a name of a second solvent.
        */
        std::string secondSolvent() const;
        
        /*!
            Sets a name of a second solvent.
        */
        //std::string setSecondSolvent(std::string);
        
        /*!
            Returns a concentration of second solvent in the component A.
        */
        double secondSolventConcentrationA() const;

        /*!
            Returns a concentration of second solvent in the component B.
        */
        double secondSolventConcentrationB() const;
        
        /*!
            Returns an interator to the first point of the gradient. The first 
            value in the pair is time and the second is concentration of
            the component B.
        */
        std::vector<std::pair<double,double> >::const_iterator beginGradient()
        const;
        
        /*!
            Returns an interator to the next to the last point of the gradient.
            The first value in the pair is time and the second is concentration
            of the component B.
        */
        std::vector<std::pair<double,double> >::const_iterator endGradient()
        const;
        
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
        gradientProfile mGradient;
        std::string mSecondSolvent;
        double mSecondSolventConcentrationA;
        double mSecondSolventConcentrationB;
};

#endif
