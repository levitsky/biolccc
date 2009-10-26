SET(theorchromo_SRCS 
   main.cpp
   chemicalbasis.cpp
   aminoacid.cpp
   terminus.cpp
   chromoconditions.cpp
   peptidemethods.cpp
)

SET(theorchromo_HDRS
   chemicalbasis.h
   aminoacid.h
   terminus.h
   chromoconditions.h
   peptidemethods.h
)

ADD_EXECUTABLE(theorchromo ${theorchromo_SRCS} ${theorchromo_HDRS})

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

SET_SOURCE_FILES_PROPERTIES(theorchromo.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(theorchromo_lib python
                theorchromo.i
                ${theorchromo_SRCS})
SWIG_LINK_LIBRARIES(theorchromo_lib ${PYTHON_LIBRARIES})

