cmake_minimum_required(VERSION 3.16)
project(lc-qt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Qt6
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

add_executable(lc-qt main.cpp)
target_link_libraries(lc-qt PRIVATE Qt6::Widgets)
