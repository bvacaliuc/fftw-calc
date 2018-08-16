cmake_minimum_required(VERSION 2.8)

# http://stackoverflow.com/questions/4029587/compiling-32-bit-vs-64-bit-project-using-cmake
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( FFTW_URL  "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4-dll64.zip" )
  set( FFTW_PATH "fftw-3.3.4-dll64.zip" )
  set( FFTW_MD5  "c840c5450f57a98578783fa633047d83" )
  set( FFTW_ARCH "x64" )
else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( FFTW_URL  "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4-dll32.zip" )
  set( FFTW_PATH "fftw-3.3.4-dll32.zip" )
  set( FFTW_MD5  "757ce6282e7cab2dfa2dd40ba46f1c1f" )
  set( FFTW_ARCH "x86" )
endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# http://stackoverflow.com/questions/8153519/how-to-automatically-download-c-dependencies-in-a-cross-platform-way-cmake
set( FFTW_HEADER "fftw3.h")
set( FFTW_LIB    "libfftw3-3.lib")
set( FFTW_DEF    "libfftw3-3.def")
get_filename_component(dir "${CMAKE_CURRENT_LIST_FILE}" PATH)

if (NOT EXISTS "${dir}/${FFTW_PATH}")
    MESSAGE( STATUS "Win32: downloading ${FFTW_URL}" )
    # https://cmake.org/pipermail/cmake/2011-December/048305.html
    file(DOWNLOAD "${FFTW_URL}" "${dir}/${FFTW_PATH}" EXPECTED_HASH MD5=${FFTW_MD5})
endif()

# http://stackoverflow.com/questions/7050997/zip-files-using-cmake
##find_package(Java)
##add_custom_command(
##    OUTPUT "${FFTW_HEADER}"
##    COMMAND
##        "${Java_JAR_EXECUTABLE}" "xfM" "${FFTW_PATH}" 
##    DEPENDS "${FFTW_PATH}"
##    )
# Need cmake 3.2+ for this
##add_custom_target(create_zip COMMAND
##    ${CMAKE_COMMAND} -E tar "xfv" "${FFTW_PATH}" --format=zip
##    )
# http://stackoverflow.com/questions/17546016/how-can-you-zip-or-unzip-from-the-command-prompt-using-only-windows-built-in-ca
##add_custom_command(
##    OUTPUT "${FFTW_HEADER}"
##    COMMAND
##        "powershell.exe" "-nologo -noprofile -command" "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('${FFTW_PATH}', '*'); }"
##    DEPENDS "${FFTW_PATH}"
##    )
# so punt, for now...
if (NOT EXISTS "${dir}/${FFTW_ARCH}/${FFTW_HEADER}")
    MESSAGE( FATAL_ERROR "Win32: please extract the contents of ${dir}/${FFTW_PATH} to ${dir}/${FFTW_ARCH} and re-run cmake" )
endif()
if (NOT EXISTS "${dir}/${FFTW_ARCH}/${FFTW_LIB}")
    MESSAGE( STATUS "Win32: generating ${FFTW_LIB} in ${dir}/${FFTW_ARCH}" )
    execute_process(
        COMMAND lib "/def:${FFTW_DEF}"
        WORKING_DIRECTORY "${dir}/${FFTW_ARCH}"
        )
endif()
