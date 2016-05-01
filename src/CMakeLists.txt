cmake_minimum_required(VERSION 2.8)

MACRO( FIND_PATH_REQUIRED arg1 arg2 )
    FIND_PATH( ${arg1} ${arg2} ${ARGN} )
    IF( "${${arg1}}" STREQUAL "${arg1}-NOTFOUND"  )
        MESSAGE( FATAL_ERROR "${arg2} not found" )
    ELSE()
        MESSAGE( STATUS "${arg2} found in ${${arg1}}" )
    ENDIF()
ENDMACRO( FIND_PATH_REQUIRED )

PROJECT( fftwcalc )
SET( fftwcalc_SRCS fftw-calc.c cpuid.c )
IF(WIN32)
    # prepare dependencies under windows
    add_subdirectory(${CMAKE_SOURCE_DIR}/../dependencies/fftw ${CMAKE_CURRENT_BINARY_DIR}/../dependencies/fftw)
ELSE()
    MESSAGE( AUTHOR_WARNING "FFTw Must already be present and installed on this system" )
ENDIF(WIN32)
# Configure to link FFTw
FIND_PATH_REQUIRED( FFTW_INC fftw3.h PATHS "${CMAKE_SOURCE_DIR}/../dependencies/fftw" PATH_SUFFIXES include )
INCLUDE_DIRECTORIES( ${FFTW_INC} )
IF(WIN32)
    MESSAGE( STATUS "Win32: using libfftw3-3.lib for linkage" )
    FIND_LIBRARY( FFTW_LIB libfftw3-3.lib REQUIRED PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
ELSEIF(UNIX)
    MESSAGE( STATUS "unix: using libfftw3.a for linkage" )
    FIND_LIBRARY( FFTW_LIB fftw3 REQUIRED PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
ELSE()
    # https://cmake.org/Wiki/CMake_Useful_Variables
    MESSAGE( WARNING "Untested for ${CMAKE_SYSTEM_NAME}:${CMAKE_CXX_COMPILER_ID}" )
    FIND_LIBRARY( FFTW_LIB fftw3 REQUIRED PATHS ${FFTW_INC} ${FFTW_INC}/../lib )
ENDIF()
# Configure to link OpenMP
FIND_PACKAGE( OpenMP )
IF(OPENMP_FOUND)
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
    SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
ENDIF(OPENMP_FOUND)
ADD_EXECUTABLE( fftw-calc ${fftwcalc_SRCS} )

# Handle linking and post-build steps
IF(WIN32)
    TARGET_LINK_LIBRARIES( fftw-calc ${LIBS} ${FFTW_LIB} )
    # http://stackoverflow.com/questions/10046114/in-cmake-how-can-i-test-if-the-compiler-is-clang
    # http://stackoverflow.com/questions/7988297/cmake-to-add-vs2010-project-custom-build-events
    IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC" OR MINGW)
        STRING( REPLACE ".lib" ".dll" FFTW_DLL ${FFTW_LIB} )
        MESSAGE( STATUS "FFTW_DLL=${FFTW_DLL}" )
        add_custom_command(TARGET fftw-calc
                           POST_BUILD
                           COMMAND ${CMAKE_COMMAND} -E copy "${FFTW_DLL}" \$\(OutDir\)
                           COMMENT "Copy ${FFTW_DLL} to OutDir")
    ENDIF()
ELSE()
    TARGET_LINK_LIBRARIES( fftw-calc ${LIBS} ${FFTW_LIB} m )
ENDIF()
