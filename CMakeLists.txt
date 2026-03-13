cmake_minimum_required(VERSION 3.26)
project(FastFuriousTransformation
        VERSION 1.0.0
        DESCRIPTION "Fast Fourier Transformation implementation"
        HOMEPAGE_URL "https://github.com/kingkybel/FastFuriousTransformation")

include(GNUInstallDirs)
include(CTest)

# Use shared settings from cmake-common submodule
include(cmake-common/cmake/DkybBuildSettings.cmake)
dkyb_apply_common_settings(CXX_STANDARD 23)

set(CMAKE_FIND_LIBRARY_PREFIXES lib)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(
        BEFORE SYSTEM ${CMAKE_SOURCE_DIR}/include
        system
)
link_directories(
        /usr/local/lib
        system
)

add_library(FastFuriousTransformation INTERFACE)
# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
target_include_directories(
        FastFuriousTransformation
        INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(FastFuriousTransformation INTERFACE cxx_std_23)

add_subdirectory(test)

set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=* -p ${CMAKE_SOURCE_DIR}/build")

file(GLOB HEADER_FILES ${CMAKE_SOURCE_DIR}/include/*.h)

# Install the header files directly under 'include/dkyb'
install(FILES ${HEADER_FILES}
        DESTINATION include/dkyb)
