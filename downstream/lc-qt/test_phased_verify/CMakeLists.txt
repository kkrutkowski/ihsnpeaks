cmake_minimum_required(VERSION 3.16)
project(phased_verify LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_AUTOMOC ON)

if(NOT DEFINED EXTERNAL_INCLUDE_DIR)
    set(EXTERNAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
endif()
if(NOT DEFINED LCQT_DIR)
    set(LCQT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
endif()
if(NOT DEFINED LCQT_BUILD_DIR)
    set(LCQT_BUILD_DIR "${LCQT_DIR}/build")
endif()

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

add_executable(phased_verify
    main_verify.cpp
    ${LCQT_DIR}/src/windows/phased_lightcurve.cpp
    ${LCQT_DIR}/src/windows/zoomed_spectrum.cpp
)

target_include_directories(phased_verify PRIVATE
    ${LCQT_DIR}/src
    ${EXTERNAL_INCLUDE_DIR}
)

# Reuse the already-compiled backend objects (loader + detrend + NuFFT)
target_link_libraries(phased_verify PRIVATE
    Qt6::Widgets
    m
    ${LCQT_BUILD_DIR}/CMakeFiles/lc_backend.dir/src/lc_periodogram.c.o
    ${LCQT_BUILD_DIR}/CMakeFiles/lc_backend.dir/upstream_src/nufft/nufft1.c.o
)
target_link_options(phased_verify PRIVATE -static -pthread)
