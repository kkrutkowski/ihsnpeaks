cmake_minimum_required(VERSION 3.16)
project(lc-qt LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_AUTOMOC ON)

# Garbage collection of unused sections to reduce binary size
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffunction-sections -fdata-sections")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffunction-sections -fdata-sections")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -Wl,--gc-sections")

# Upstream source/include locations (overridable; build.sh mounts them into the container)
if(NOT DEFINED EXTERNAL_INCLUDE_DIR)
    set(EXTERNAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
endif()
if(NOT DEFINED UPSTREAM_SRC_DIR)
    set(UPSTREAM_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
endif()

# Microarchitecture for the backend (math + NuFFT + scaling constants). Empty = toolchain
# default ISA. Dispatch-ready: a future multi-variant build compiles lc_backend once per
# value of LC_BACKEND_MARCH with per-variant symbol suffixes + a CPUID dispatcher.
set(LC_BACKEND_MARCH "" CACHE STRING "Backend microarchitecture (e.g. x86-64, x86-64-v3); empty = toolchain default")

# ---------------------------------------------------------------------------------------
# Generate scaling.h at build time (mirrors the ihsnpeaks scaling_gen step) instead of
# shipping the generic constants. scaling.c + nufft1.c are compiled into a small static
# helper that is run to emit tuned F_ALPHA / F_PSWF*_BETA / F_PSWF*_GAMMA into the build dir.
# ---------------------------------------------------------------------------------------
set(SCALING_NUFFT_OBJ "${CMAKE_BINARY_DIR}/nufft1.scaling.o")
set(SCALING_GEN "${CMAKE_BINARY_DIR}/scaling_gen")
set(SCALING_HEADER "${CMAKE_BINARY_DIR}/scaling.h")

set(SCALING_FLAGS -std=gnu23 -D_GNU_SOURCE -DHAS_MIMALLOC=0 -DMAX_TWIDDLE_REUSE=8 -O3 -ffast-math -fno-sanitize=all)
if(LC_BACKEND_MARCH)
    list(APPEND SCALING_FLAGS -march=${LC_BACKEND_MARCH} -mtune=generic)
endif()

add_custom_command(
    OUTPUT "${SCALING_HEADER}"
    COMMAND ${CMAKE_C_COMPILER} ${SCALING_FLAGS} -I${UPSTREAM_SRC_DIR}/nufft
            -c ${UPSTREAM_SRC_DIR}/nufft/nufft1.c -o ${SCALING_NUFFT_OBJ}
    COMMAND ${CMAKE_C_COMPILER} ${SCALING_FLAGS} -I${UPSTREAM_SRC_DIR}/nufft
            ${UPSTREAM_SRC_DIR}/nufft/scaling.c ${SCALING_NUFFT_OBJ}
            -static -lm -o ${SCALING_GEN}
    COMMAND ${SCALING_GEN} ${SCALING_HEADER}
    DEPENDS ${UPSTREAM_SRC_DIR}/nufft/nufft1.c ${UPSTREAM_SRC_DIR}/nufft/scaling.c ${UPSTREAM_SRC_DIR}/nufft/nufft1.h
    COMMENT "Generating scaling.h via scaling_gen"
    VERBATIM
)
add_custom_target(scaling_header DEPENDS "${SCALING_HEADER}")

# ---------------------------------------------------------------------------------------
# Backend OBJECT library: a single bridge object (lc_periodogram.o) plus a separate NuFFT
# object (nufft1.o). The LC_API markers in the headers document the intended public entry
# points for a future multi-microarchitecture dispatch (mirrors the ihsnpeaks release model);
# no -fvisibility=hidden is used, so all symbols are exported and link normally.
# ---------------------------------------------------------------------------------------
add_library(lc_backend OBJECT
    src/lc_periodogram.c
    ${UPSTREAM_SRC_DIR}/nufft/nufft1.c
)
add_dependencies(lc_backend scaling_header)

target_compile_definitions(lc_backend PRIVATE _GNU_SOURCE HAS_MIMALLOC=0 MAX_TWIDDLE_REUSE=8)
target_compile_options(lc_backend PRIVATE -fno-sanitize=all -pthread)
if(LC_BACKEND_MARCH)
    target_compile_options(lc_backend PRIVATE -march=${LC_BACKEND_MARCH} -mtune=generic)
endif()
target_include_directories(lc_backend PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EXTERNAL_INCLUDE_DIR}
    ${UPSTREAM_SRC_DIR}/nufft
    ${UPSTREAM_SRC_DIR}
    ${CMAKE_BINARY_DIR}
)

# Find Qt6
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

add_executable(lc-qt
    src/main.cpp
    src/windows/customize_labels.cpp
    src/windows/classification_stats.cpp
    src/windows/lightcurve_plot.cpp
    src/windows/spectrum_plot.cpp
    src/windows/customize_period_search.cpp
)

# Global include paths (safe for moc)
target_include_directories(lc-qt PRIVATE
    src
    ${EXTERNAL_INCLUDE_DIR}
)

target_link_libraries(lc-qt PRIVATE Qt6::Widgets m lc_backend)
target_link_options(lc-qt PRIVATE -static -pthread)
