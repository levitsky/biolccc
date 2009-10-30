// pyBioLCCC.i - SWIG interface
%module pyBioLCCC 
%include "std_string.i"

%{
#include "aminoacid.h"
#include "terminus.h"
#include "chemicalbasis.h"
#include "gradientpoint.h"
#include "gradient.h"
#include "chromoconditions.h"
#include "BioLCCC.h"
%}

// Parse the original header file
%include "aminoacid.h"
%include "terminus.h"
%include "chemicalbasis.h"
%include "gradientpoint.h"
%include "gradient.h"
%include "chromoconditions.h"
%include "BioLCCC.h"

// Instantiate some templates

