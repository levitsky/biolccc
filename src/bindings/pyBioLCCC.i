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
%template(StringChemicalGroupMap) std::map<std::string,BioLCCC::ChemicalGroup>;

%extend BioLCCC::ChemicalGroup{
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

            output_dict['chemicalGroups'] = {}
            for label, chemical_group in self.chemicalGroups().items():
                output_dict['chemicalGroups'][label] = (
                    chemical_group.__getstate__())

            return output_dict

        def __setstate__(self, chembasis_dict):
            self.setModel(chembasis_dict['model'])
            self.setSegmentLength(chembasis_dict['segmentLength'])
            self.setPersistentLength(chembasis_dict['persistentLength'])
            self.setAdsorbtionLayerWidth(chembasis_dict['adsorbtionLayerWidth'])
            self.setSecondSolventBindEnergy(
                chembasis_dict['secondSolventBindEnergy'])
            self.clearChemicalGroups()
            for chemical_group_dict in chembasis_dict['chemicalGroups'].values():
                chemical_group = ChemicalGroup(
                    chemical_group_dict['name'],
                    chemical_group_dict['label'],
                    chemical_group_dict['bindEnergy'],
                    chemical_group_dict['averageMass'],
                    chemical_group_dict['monoisotopicMass'])
                self.addChemicalGroup(chemical_group)

        def min_inf(self):
            output_dict = {}
            output_dict['model'] = self.model()
            output_dict['segmentLength'] = self.segmentLength()
            output_dict['persistentLength'] = self.persistentLength()
            output_dict['adsorbtionLayerWidth'] = self.adsorbtionLayerWidth()
            output_dict['secondSolventBindEnergy'] = \
                self.secondSolventBindEnergy()

            for label, chemical_group in self.chemicalGroups().items():
                output_dict[label] = chemical_group.bindEnergy()
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
                if key in self.chemicalGroups():
                    self.setChemicalGroupBindEnergy(key, value) 
                else:
                    self.addChemicalGroup(ChemicalGroup("", key, value))

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
%include "chemicalgroup.h"
%include "chemicalbasis.h"
%include "gradientpoint.h"
%include "gradient.h"
%include "chromoconditions.h"
%include "BioLCCC.h"

// Instantiate some templates

