// pyBioLCCC.i - SWIG interface
%module pyBioLCCC 

%feature("autodoc", "0");

%{
#include "auxiliary.hpp"
#include "chemicalgroup.h"
#include "chemicalbasis.h"
#include "gradientpoint.h"
#include "gradient.h"
#include "chromoconditions.h"
#include "BioLCCC.h"
%}

%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"
%template(GradientPointVector) std::vector<BioLCCC::GradientPoint>;
%template(StringAminoacidMap) std::map<std::string,BioLCCC::ChemicalGroup>;

%extend BioLCCC::Aminoacid {
    %insert("python") %{
        def __getstate__(self):
            output_dict = {}
            output_dict['name'] = self.name()
            output_dict['label'] = self.label()
            output_dict['bindEnergy'] = self.bindEnergy()
            output_dict['averageMass'] = self.averageMass()
            output_dict['monoisotopicMass'] = self.monoisotopicMass()
            return output_dict
    %}
}

%extend BioLCCC::Terminus {
    %insert("python") %{
        def __getstate__(self):
            output_dict = {}
            output_dict['name'] = self.name()
            output_dict['label'] = self.label()
            output_dict['bindEnergy'] = self.bindEnergy()
            output_dict['averageMass'] = self.averageMass()
            output_dict['monoisotopicMass'] = self.monoisotopicMass()
            return output_dict
    %}
}

%extend BioLCCC::ChemicalBasis {
    %insert("python") %{
        def __str__(self):
            return str(self.min_inf()).replace(',', ',\n')

        def __getstate__(self):
            output_dict = {}
            output_dict['model'] = self.model()
            output_dict['segmentLength'] = self.segmentLength()
            output_dict['persistentLength'] = self.persistentLength()
            output_dict['adsorbtionLayerWidth'] = self.adsorbtionLayerWidth()
            output_dict['secondSolventBindEnergy'] = \
                self.secondSolventBindEnergy()

            output_dict['aminoacids'] = {}
            for label, aminoacid in self.aminoacids().items():
                output_dict['aminoacids'][label] = aminoacid.__getstate__()

            output_dict['CTermini'] = {}
            for label, CTerminus in self.CTermini().items():
                output_dict['CTermini'][label] = CTerminus.__getstate__()

            output_dict['NTermini'] = {}
            for label, NTerminus in self.NTermini().items():
                output_dict['NTermini'][label] = NTerminus.__getstate__()
            return output_dict

        def __setstate__(self, chembasis_dict):
            self.setModel(chembasis_dict['model'])
            self.setSegmentLength(chembasis_dict['segmentLength'])
            self.setPersistentLength(chembasis_dict['persistentLength'])
            self.setAdsorbtionLayerWidth(chembasis_dict['adsorbtionLayerWidth'])
            self.setSecondSolventBindEnergy(
                chembasis_dict['secondSolventBindEnergy'])
            self.clearAminoacids()
            for aminoacid_dict in chembasis_dict['aminoacids'].values():
                aminoacid = Aminoacid(aminoacid_dict['name'],
                                      aminoacid_dict['label'],
                                      aminoacid_dict['bindEnergy'],
                                      aminoacid_dict['averageMass'],
                                      aminoacid_dict['monoisotopicMass'])
                self.addAminoacid(aminoacid)
            self.clearCTermini()
            for CTerminus_dict in chembasis_dict['CTermini'].values():
                CTerminus = Terminus(CTerminus_dict['name'],
                                     CTerminus_dict['label'],
                                     CTerminus_dict['bindEnergy'],
                                     CTerminus_dict['averageMass'],
                                     CTerminus_dict['monoisotopicMass'])
                self.addCTerminus(CTerminus)
            self.clearNTermini()
            for NTerminus_dict in chembasis_dict['NTermini'].values():
                NTerminus = Terminus(NTerminus_dict['name'],
                                     NTerminus_dict['label'],
                                     NTerminus_dict['bindEnergy'],
                                     NTerminus_dict['averageMass'],
                                     NTerminus_dict['monoisotopicMass'])
                self.addNTerminus(NTerminus)

        def min_inf(self):
            output_dict = {}
            output_dict['model'] = self.model()
            output_dict['segmentLength'] = self.segmentLength()
            output_dict['persistentLength'] = self.persistentLength()
            output_dict['adsorbtionLayerWidth'] = self.adsorbtionLayerWidth()
            output_dict['secondSolventBindEnergy'] = \
                self.secondSolventBindEnergy()

            for label, aminoacid in self.aminoacids().items():
                output_dict[label] = aminoacid.bindEnergy()

            for label, CTerminus in self.CTermini().items():
                output_dict[label] = CTerminus.bindEnergy()

            for label, NTerminus in self.NTermini().items():
                output_dict[label] = NTerminus.bindEnergy()
            return output_dict

        def set_min_inf_element(self, key, value):
            if key == 'model':
                self.setModel(value)
            elif key == 'segmentLength':
                self.setSegmentLength(value)
            elif key == 'persistentLength':
                self.setPersistentLength(value)
            elif key == 'adsorbtionLayerWidth':
                self.setAdsorbtionLayerWidth(value)
            elif key == 'secondSolventBindEnergy':
                self.setSecondSolventBindEnergy(value)
            else:
                if key.endswith('-'):
                    if key in self.NTermini():
                        self.setNTerminusBindEnergy(key, value) 
                    else:
                        self.addNTerminus(Terminus("",
                                                   key, 
                                                   value))
                elif key.startswith('-'):
                    if key in self.CTermini():
                        self.setCTerminusBindEnergy(key, value) 
                    else:
                        self.addCTerminus(Terminus("",
                                                   key, 
                                                   value))
                else:
                    if key in self.aminoacids():
                        self.setAminoacidBindEnergy(key, value) 
                    else:
                        self.addAminoacid(Aminoacid("",
                                                    key, 
                                                    value))

        def set_min_inf(self, min_inf_dict):
            for key, value in min_inf_dict.items():
                self.set_min_inf_element(key, value)
    %}
};

%extend BioLCCC::GradientPoint{
    %insert("python") %{
        def __str__(self):
            return '(%f,%f)' % (self.time(), self.concentrationB())
    %}
};

%extend BioLCCC::ChromoConditions{
    %insert("python") %{
        def __str__(self):
            return str(self.__getstate__())

        def __getstate__(self):
            output_dict = {}
            output_dict['columnLength'] = self.columnLength()
            output_dict['columnDiameter'] = self.columnDiameter()
            output_dict['columnPoreSize'] = self.columnPoreSize()
            output_dict['gradient'] = [eval(str(i)) for i in self.gradient()]
            output_dict['secondSolventConcentrationA'] = \
                self.secondSolventConcentrationA()
            output_dict['secondSolventConcentrationB'] = \
                self.secondSolventConcentrationB()
            output_dict['delayTime'] = self.delayTime()
            output_dict['flowRate'] = self.flowRate()
            output_dict['dV'] = self.dV()
            output_dict['calibrationParameter'] = self.calibrationParameter()
            output_dict['columnVpToVtot'] = self.columnVpToVtot()
            output_dict['columnPorosity'] = self.columnPorosity()
            output_dict['temperature'] = self.temperature()
            return output_dict

        def __setstate__(self, input_dict):
            self.setColumnLength(input_dict['columnLength'])
            self.setColumnDiameter(input_dict['columnDiameter'])
            self.setColumnPoreSize(input_dict['columnPoreSize'])

            input_gradient = self.gradient();
            input_gradient.clear();
            for (time, concentrationB) in input_dict['gradient']:
                input_gradient.addPoint(time, concentrationB)
            self.setGradient(input_gradient)
            self.setSecondSolventConcentrationA(
                input_dict['secondSolventConcentrationA'])
            self.setSecondSolventConcentrationB(
                input_dict['secondSolventConcentrationB'])
            self.setDelayTime(input_dict['delayTime'])
            self.setFlowRate(input_dict['flowRate'])
            self.setDV(input_dict['dV'])
            self.setCalibrationParameter(input_dict['calibrationParameter'])
            self.setColumnVpToVtot(input_dict['columnVpToVtot'])
            self.setColumnPorosity(input_dict['columnPorosity'])
            self.setTemperature(input_dict['temperature'])
    %}
};

// Parse the original header file
%include "auxiliary.hpp"
%include "aminoacid.h"
%include "terminus.h"
%include "chemicalbasis.h"
%include "gradientpoint.h"
%include "gradient.h"
%include "chromoconditions.h"
%include "BioLCCC.h"

// Instantiate some templates

