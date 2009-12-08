INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(BioLCCC_SRCS 
   aminoacid.cpp
   terminus.cpp
   chemicalbasis.cpp
   gradientpoint.cpp
   gradient.cpp
   chromoconditions.cpp
   BioLCCC.cpp
)

SET(BioLCCC_HDRS
   auxiliary.hpp
   aminoacid.h
   terminus.h
   chemicalbasis.h
   gradientpoint.h
   gradient.h
   chromoconditions.h
   BioLCCC.h
)

ADD_EXECUTABLE(BioLCCC_bin main.cpp ${BioLCCC_SRCS} ${BioLCCC_HDRS})
ADD_EXECUTABLE(BioLCCC_test tests.cpp ${BioLCCC_SRCS} tests.h ${BioLCCC_HDRS})

IF(UNIX)
    add_definitions(-Wall -O2)
ENDIF(UNIX)

FIND_PACKAGE(SWIG)
IF(SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
    SET(CMAKE_SWIG_FLAGS "")

    FIND_PACKAGE(PythonLibs)
    IF(PythonLibs_FOUND)
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

        SET_SOURCE_FILES_PROPERTIES(pyBioLCCC.i PROPERTIES CPLUSPLUS ON)
        SWIG_ADD_MODULE(pyBioLCCC python
                        pyBioLCCC.i
                        ${BioLCCC_SRCS})
        SWIG_LINK_LIBRARIES(pyBioLCCC ${PYTHON_LIBRARIES})
    ENDIF(PythonLibs_FOUND)

    SET_SOURCE_FILES_PROPERTIES(csBioLCCC.i PROPERTIES CPLUSPLUS ON)
    SWIG_ADD_MODULE(csBioLCCC csharp 
                    csBioLCCC.i
                    ${BioLCCC_SRCS})
ENDIF(SWIG_FOUND)
#SWIG_LINK_LIBRARIES(csBioLCCC ${PYTHON_LIBRARIES})
