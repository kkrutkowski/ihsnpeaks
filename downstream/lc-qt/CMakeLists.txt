cmake_minimum_required(VERSION 3.16)
project(lc-qt LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 23)
set(CMAKE_AUTOMOC ON)

# Garbage collection of unused sections to reduce binary size
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffunction-sections -fdata-sections")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffunction-sections -fdata-sections")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -Wl,--gc-sections")


# Find Qt6
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

add_executable(lc-qt
    src/main.cpp
    src/lc_readout.c
    src/windows/customize_labels.cpp
    src/windows/classification_stats.cpp
    src/windows/lightcurve_plot.cpp
)

if(NOT DEFINED EXTERNAL_INCLUDE_DIR)
    set(EXTERNAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
endif()

if(NOT DEFINED UPSTREAM_SRC_DIR)
    set(UPSTREAM_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
endif()

# Global include paths (safe for moc)
target_include_directories(lc-qt PRIVATE
    src
    ${EXTERNAL_INCLUDE_DIR}
)

# Upstream include paths ONLY for the C bridge (avoids moc header conflicts)
# Note: readout.h uses relative includes, so we only need nufft and src root
set_source_files_properties(src/lc_readout.c PROPERTIES COMPILE_FLAGS
    "-std=gnu23 -I${UPSTREAM_SRC_DIR}/nufft -I${UPSTREAM_SRC_DIR}")

target_link_libraries(lc-qt PRIVATE Qt6::Widgets m)
target_link_options(lc-qt PRIVATE -static)
