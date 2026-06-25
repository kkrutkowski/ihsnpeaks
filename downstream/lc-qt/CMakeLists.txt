cmake_minimum_required(VERSION 3.16)
project(lc-qt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

# Garbage collection of unused sections to reduce binary size
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffunction-sections -fdata-sections")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -Wl,--gc-sections")


# Find Qt6
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

add_executable(lc-qt
    src/main.cpp
    src/windows/customize_labels.cpp
    src/windows/classification_stats.cpp
)

if(NOT DEFINED EXTERNAL_INCLUDE_DIR)
    set(EXTERNAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
endif()

target_include_directories(lc-qt PRIVATE src ${EXTERNAL_INCLUDE_DIR})
target_link_libraries(lc-qt PRIVATE Qt6::Widgets)
target_link_options(lc-qt PRIVATE -static)
