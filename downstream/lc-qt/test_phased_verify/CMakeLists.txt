cmake_minimum_required(VERSION 3.16)
project(phased_verify LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_AUTOMOC ON)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

if(NOT DEFINED EXTERNAL_INCLUDE_DIR)
    set(EXTERNAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
endif()
if(NOT DEFINED UPSTREAM_SRC_DIR)
    set(UPSTREAM_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../src")
endif()
if(NOT DEFINED LCQT_DIR)
    set(LCQT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
endif()

# Generate scaling.h
set(SCALING_NUFFT_OBJ "${CMAKE_BINARY_DIR}/nufft1.scaling.o")
set(SCALING_GEN "${CMAKE_BINARY_DIR}/scaling_gen")
set(SCALING_HEADER "${CMAKE_BINARY_DIR}/scaling.h")
set(SCALING_FLAGS -std=gnu23 -D_GNU_SOURCE -DHAS_MIMALLOC=0 -DMAX_TWIDDLE_REUSE=8 -O3 -ffast-math -fno-sanitize=all)

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

add_library(lc_backend OBJECT
    ${LCQT_DIR}/src/lc_period.c
    ${UPSTREAM_SRC_DIR}/nufft/nufft1.c
)
add_dependencies(lc_backend scaling_header)
target_compile_definitions(lc_backend PRIVATE _GNU_SOURCE HAS_MIMALLOC=0 MAX_TWIDDLE_REUSE=8)
target_compile_options(lc_backend PRIVATE -fno-sanitize=all -pthread)
target_include_directories(lc_backend PRIVATE
    ${LCQT_DIR}/src
    ${EXTERNAL_INCLUDE_DIR}
    ${UPSTREAM_SRC_DIR}/nufft
    ${UPSTREAM_SRC_DIR}
    ${CMAKE_BINARY_DIR}
)

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

target_link_libraries(phased_verify PRIVATE
    Qt6::Widgets
    m
    lc_backend
)
target_link_options(phased_verify PRIVATE -static -pthread)

add_executable(test_render_png
    test_render_png.cpp
    ${LCQT_DIR}/src/windows/phased_lightcurve.cpp
)
target_include_directories(test_render_png PRIVATE
    ${LCQT_DIR}/src
    ${EXTERNAL_INCLUDE_DIR}
)
target_link_libraries(test_render_png PRIVATE
    Qt6::Widgets
    m
    lc_backend
)
target_link_options(test_render_png PRIVATE -static -pthread)
