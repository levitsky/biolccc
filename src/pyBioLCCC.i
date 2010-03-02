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

