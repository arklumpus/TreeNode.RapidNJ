# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)

project ("rapidNJWrapper")

# Include sub-projects.
add_subdirectory ("rapidNJWrapper")
