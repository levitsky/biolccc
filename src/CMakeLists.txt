SET(BioLCCC_SRCS 
   main.cpp
   chemicalbasis.cpp
   aminoacid.cpp
   terminus.cpp
   chromoconditions.cpp
   peptidemethods.cpp
)

SET(BioLCCC_HDRS
   chemicalbasis.h
   aminoacid.h
   terminus.h
   chromoconditions.h
   peptidemethods.h
)

ADD_EXECUTABLE(BioLCCC_bin ${BioLCCC_SRCS} ${BioLCCC_HDRS})

IF(UNIX)
    add_definitions(-Wall -O2)
ENDIF(UNIX)

# This example shows how to use python
# Currently these languages have been tested:
#  perl tcl ruby php4 pike

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(pyBioLCCC.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(pyBioLCCC python
                pyBioLCCC.i
                ${BioLCCC_SRCS})
SWIG_LINK_LIBRARIES(pyBioLCCC ${PYTHON_LIBRARIES})

