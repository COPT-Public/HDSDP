cmake_minimum_required(VERSION 3.21)
# set(CMAKE_OSX_ARCHITECTURES x86_64)
project(hdsdp LANGUAGES C)
include(GNUInstallDirs)
include(CTest)

#-------------------------------
# Options and configurations
#-------------------------------
option(LINSYS_PARDISO "Enable MKL Pardiso support (requires MKL)" OFF)

# Add compile definitions based on user input
if(LINSYS_PARDISO)
    add_compile_definitions(LINSYS_PARDISO)
endif()

add_compile_definitions(HEADERPATH)
# ddot_
# add_compile_definitions(UNDERBLAS)    
# DDOT
# add_compile_definitions(CAPBLAS)
# DDOT_
# add_compile_definitions(CAPUNDERBLAS)

#-------------------------------
# Find BLAS/LAPACK
#-------------------------------
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(FATAL_ERROR "BLAS and LAPACK are required but not found.")
endif()

#-------------------------------
# If LINSYS_PARDISO is enabled, link MKL
# Adjust this section for your MKL installation.
#-------------------------------
if(LINSYS_PARDISO)
    # Attempt to find MKL via environment variables or specify manually
    # You may need something like:
    # set(MKL_ROOT "/path/to/mkl")
    # include_directories("${MKL_ROOT}/include")
    # list(APPEND MKL_LIBS "${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so" "${MKL_ROOT}/lib/intel64/libmkl_core.so" "${MKL_ROOT}/lib/intel64/libmkl_sequential.so")

    # The exact MKL libraries to link depend on your system and MKL version.
    # For demonstration, we'll assume user sets MKL_LIBS externally or knows how to modify this part.
    if(WIN32)
        # Example for Windows (adjust .lib files as needed):
        # list(APPEND MKL_LIBS "mkl_intel_lp64.lib" "mkl_core.lib" "mkl_sequential.lib")
    elseif(APPLE)
        # macOS Example:
        list(APPEND MKL_LIBS "/opt/intel/oneapi/mkl/2021.2.0/lib/libmkl_intel_lp64.dylib" "/opt/intel/oneapi/mkl/2021.2.0/lib/libmkl_core.dylib" "/opt/intel/oneapi/mkl/2021.2.0/lib/libmkl_sequential.dylib")
    else()
        # Linux Example:
        # list(APPEND MKL_LIBS "/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_intel_lp64.so" "/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_core.so" "/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_sequential.so")
    endif()
endif()

#-------------------------------
# Gather sources for libraries
#
# We want all sources outside test/, excluding memwatch.*
# The directory structure is:
#  - external
#  - interface
#  - linalg
#
# Also exclude main.c (it has a main function or might cause conflicts), and exclude memwatch.
#-------------------------------

file(GLOB EXTERNAL_SOURCES
    external/*.c
)
file(GLOB INTERFACE_SOURCES
    interface/*.c
)
file(GLOB LINALG_SOURCES
    linalg/*.c
)

# Exclude memwatch from external sources
list(FILTER EXTERNAL_SOURCES EXCLUDE REGEX "memwatch\\.c")

# Exclude main.c if present in top-level (to avoid multiple mains or unwanted code in the library)
set(TOP_LEVEL_SOURCES "")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/main.c")
    set(TOP_LEVEL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.c")
endif()

# Exclude main.c from library sources (if it exists)
list(FILTER TOP_LEVEL_SOURCES EXCLUDE REGEX "main\\.c")

set(LIB_SOURCES
    ${EXTERNAL_SOURCES}
    ${INTERFACE_SOURCES}
    ${LINALG_SOURCES}
)

#-------------------------------
# Create Shared and Static libraries named hdsdp
#
# We'll create two targets: hdsdp_shared and hdsdp_static
# Then set their OUTPUT_NAME to "hdsdp"
#-------------------------------
add_library(hdsdp_shared SHARED ${LIB_SOURCES})
target_include_directories(hdsdp_shared PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hdsdp>
)
set_target_properties(hdsdp_shared PROPERTIES OUTPUT_NAME "hdsdp")
target_link_libraries(hdsdp_shared PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(LINSYS_PARDISO AND MKL_LIBS)
    target_link_libraries(hdsdp_shared PRIVATE ${MKL_LIBS})
endif()

add_library(hdsdp_static STATIC ${LIB_SOURCES})
target_include_directories(hdsdp_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hdsdp>
)
set_target_properties(hdsdp_static PROPERTIES OUTPUT_NAME "hdsdp")
target_link_libraries(hdsdp_static PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(LINSYS_PARDISO AND MKL_LIBS)
    target_link_libraries(hdsdp_static PRIVATE ${MKL_LIBS})
endif()

#-------------------------------
# Create the executable sdpasolve from tests/sdpasolve.c
#
# The sdpasolve binary presumably uses the hdsdp library.
# We explicitly include sdpasolve.c from the tests directory.
#-------------------------------
if(NOT WIN32)
    add_executable(sdpasolve tests/sdpasolve.c)
    target_include_directories(sdpasolve PUBLIC
        ${CMAKE_SOURCE_DIR}
    )
    target_link_libraries(sdpasolve PRIVATE hdsdp_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} m)
    if(LINSYS_PARDISO AND MKL_LIBS)
        target_link_libraries(sdpasolve PRIVATE ${MKL_LIBS})
    endif()

    # If you need C99 standard:
    set_property(TARGET hdsdp_shared hdsdp_static sdpasolve PROPERTY C_STANDARD 99)
endif()

install(TARGETS hdsdp_shared hdsdp_static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
if(TARGET sdpasolve)
    install(TARGETS sdpasolve RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# The headers refer to one another by interface/, linalg/, and external/ paths.
install(DIRECTORY interface linalg external
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdsdp
    FILES_MATCHING PATTERN "*.h"
)

if(BUILD_TESTING AND TARGET sdpasolve)
    add_test(NAME hdsdp_truss1
        COMMAND ${CMAKE_COMMAND}
            "-DSOLVER=$<TARGET_FILE:sdpasolve>"
            "-DINSTANCE=${CMAKE_CURRENT_SOURCE_DIR}/examples/truss1.dat-s"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_sdpasolve.cmake"
    )
endif()
