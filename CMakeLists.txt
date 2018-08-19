########################################################################
# Project setup
########################################################################
cmake_minimum_required(VERSION 2.8)
PROJECT( fftwcalc )

MACRO( FIND_PATH_REQUIRED arg1 arg2 )
    FIND_PATH( ${arg1} ${arg2} ${ARGN} )
    IF( "${${arg1}}" STREQUAL "${arg1}-NOTFOUND"  )
        MESSAGE( FATAL_ERROR "${arg2} not found" )
    ELSE()
        MESSAGE( STATUS "${arg2} found in ${${arg1}}" )
    ENDIF()
ENDMACRO( FIND_PATH_REQUIRED )

########################################################################
# Dependencies
########################################################################

IF(WIN32)
    # prepare dependencies under windows
    add_subdirectory(${PROJECT_SOURCE_DIR}/dependencies/fftw ${PROJECT_SOURCE_DIR}/dependencies/fftw)
ELSE()
    MESSAGE( AUTHOR_WARNING "FFTw Must already be present and installed on this system" )
ENDIF(WIN32)

# Configure to link FFTw
FIND_PATH_REQUIRED( FFTW_INC fftw3.h PATHS "${PROJECT_SOURCE_DIR}/dependencies/fftw/x64" "${PROJECT_SOURCE_DIR}/dependencies/fftw/x86" PATH_SUFFIXES include )
INCLUDE_DIRECTORIES( ${FFTW_INC} )
IF(WIN32)
    MESSAGE( STATUS "Win32: using libfftw3-3.lib for linkage" )
    FIND_LIBRARY( FFTW_LIB libfftw3-3.lib REQUIRED PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
    STRING( REPLACE ".lib" ".dll" FFTW_DLL ${FFTW_LIB} )
ELSEIF(UNIX)
    MESSAGE( STATUS "unix: using libfftw3.a for linkage" )
    FIND_LIBRARY( FFTW_LIB fftw3 REQUIRED PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
    MESSAGE( STATUS "unix: ${FFTW_LIB}" )
ELSE()
    # https://cmake.org/Wiki/CMake_Useful_Variables
    MESSAGE( WARNING "Untested for ${CMAKE_SYSTEM_NAME}:${CMAKE_CXX_COMPILER_ID}" )
    FIND_LIBRARY( FFTW_LIB fftw3 REQUIRED PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
    MESSAGE( STATUS "unix: ${FFTW_LIB}" )
ENDIF()
list(APPEND LIBS ${FFTW_LIB})

#link with libm when available
FIND_LIBRARY(
    MATH_LIBRARIES NAMES m
    PATHS /usr/lib /usr/lib64
)
IF(MATH_LIBRARIES)
    list(APPEND LIBS ${MATH_LIBRARIES})
ENDIF()

# Configure to link OpenMP
FIND_PACKAGE( OpenMP )
IF(UNIX)
    MESSAGE( STATUS "unix: checking for libfftw3_omp.a" )
    FIND_LIBRARY( FFTW_LIB_OMP fftw3_omp OPTIONAL PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
    list(APPEND LIBS ${FFTW_LIB_OMP})
ENDIF()
IF(OPENMP_FOUND AND EXISTS ${FFTW_LIB_OMP})
    MESSAGE( STATUS "openmp: found and using ${FFTW_LIB_OMP}" )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
    SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
ELSE()
    MESSAGE( STATUS "openmp: found but libfftw3_omp.a is not, so IGNORING it" )
ENDIF()

########################################################################
# Build executables
########################################################################

ADD_EXECUTABLE( fftw-calc src/fftw-calc.c src/cpuid.c )
TARGET_LINK_LIBRARIES( fftw-calc ${LIBS} )

ADD_EXECUTABLE( fftw-io src/fftw-io.c )
TARGET_LINK_LIBRARIES( fftw-io ${LIBS} )

# Handle linking and post-build steps
IF(WIN32)
    # http://stackoverflow.com/questions/10046114/in-cmake-how-can-i-test-if-the-compiler-is-clang
    # http://stackoverflow.com/questions/7988297/cmake-to-add-vs2010-project-custom-build-events
    IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC" OR MINGW)
        MESSAGE( STATUS "FFTW_DLL=${FFTW_DLL}" )
        add_custom_command(TARGET fftw-calc
                           POST_BUILD
                           COMMAND ${CMAKE_COMMAND} -E copy "${FFTW_DLL}" \$\(OutDir\)
                           COMMENT "Copy ${FFTW_DLL} to OutDir")
    ENDIF()
ENDIF()
