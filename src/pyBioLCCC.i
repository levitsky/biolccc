// theorchromo.i - SWIG interface
%module pyBioLCCC 
%include "std_string.i"

%{
#include "aminoacid.h"
#include "terminus.h"
#include "chromoconditions.h"
#include "chemicalbasis.h"
#include "peptidemethods.h"
%}

// Parse the original header file
%include "aminoacid.h"
%include "terminus.h"
%include "chromoconditions.h"
%include "chemicalbasis.h"
%include "peptidemethods.h"

// Instantiate some templates

