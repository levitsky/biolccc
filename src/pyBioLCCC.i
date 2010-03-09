// pyBioLCCC.i - SWIG interface
%module pyBioLCCC 

%{
#include "auxiliary.hpp"
#include "aminoacid.h"
#include "terminus.h"
#include "chemicalbasis.h"
#include "gradientpoint.h"
#include "gradient.h"
#include "chromoconditions.h"
#include "BioLCCC.h"
%}

%include "std_string.i"
%include "std_map.i"

%extend BioLCCC::ChemicalBasis {
    char *__str__() {
        static char tmp[1024];
        sprintf(tmp, "{'model': '%s'", $self->model().c_str());
        sprintf(tmp, "%s,\n'segmentLength': %.6g", 
                tmp,$self->segmentLength());
        sprintf(tmp, "%s,\n'persistentLength': %d", 
                tmp,$self->persistentLength());
        sprintf(tmp, "%s,\n'adsorbtionLayerWidth': %.6g",
                tmp, $self->adsorbtionLayerWidth());
        sprintf(tmp, "%s,\n'secondSolventBindEnergy': %.6g",
                tmp, $self->secondSolventBindEnergy());
        for (std::map<std::string, BioLCCC::Aminoacid>::const_iterator i =
                $self->aminoacids().begin();
             i != $self->aminoacids().end();
             i++) {
            sprintf(tmp, "%s,\n'%s': %.6g",
                    tmp, i->second.label().c_str(),i->second.bindEnergy());
        }

        for (std::map<std::string, BioLCCC::Terminus>::const_iterator i =
                $self->NTermini().begin();
             i != $self->NTermini().end();
             i++) {
            sprintf(tmp, "%s,\n'%s': %.6g",
                    tmp, i->second.label().c_str(),i->second.bindEnergy());
        }

        for (std::map<std::string, BioLCCC::Terminus>::const_iterator i =
                $self->CTermini().begin();
             i != $self->CTermini().end();
             i++) {
            sprintf(tmp, "%s,\n'%s': %.6g",
                    tmp, i->second.label().c_str(),i->second.bindEnergy());
        }

        sprintf(tmp, "%s}",tmp);

        return tmp;
    }
};

%pythoncode %{
    def chemicalBasis_from_dict(chembasis_dict):
        chembasis = ChemicalBasis()
        for key, val in chembasis_dict.items():
            if key == "model":
                chembasis.setModel(val)
            elif key == "segmentLength":
                chembasis.setSegmentLength(val)
            elif key == "persistentLength":
                chembasis.setPersistentLength(val)
            elif key == "adsorbtionLayerWidth":
                chembasis.setAdsorbtionLayerWidth(val)
            elif key == "secondSolventBindEnergy":
                chembasis.setSecondSolventBindEnergy(val)
            elif key.endswith("-"):
                chembasis.setNTerminusBindEnergy(key, val)
            elif key.startswith("-"):
                chembasis.setCTerminusBindEnergy(key, val)
            else:
                chembasis.setAminoacidBindEnergy(key, val)
        return chembasis
%}

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

