cmake_minimum_required(VERSION 3.15)
project(crm_bindings)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find packages
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

# Include your existing CRM library
# ADJUST THESE PATHS TO YOUR ACTUAL CRM LOCATION
set(CRM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../../CRM/include" CACHE PATH "CRM include directory")
set(CRM_LIB_DIR "${CMAKE_SOURCE_DIR}/../../../CRM/build" CACHE PATH "CRM library directory")

# Python module
pybind11_add_module(crm_cpp crm_bindings.cpp)

target_include_directories(crm_cpp PRIVATE 
    ${CRM_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_directories(crm_cpp PRIVATE ${CRM_LIB_DIR})

# Link your CRM libraries
target_link_libraries(crm_cpp PRIVATE CRMCPPLib)

# Optimization flags
target_compile_options(crm_cpp PRIVATE -O3 -march=native)

# OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(crm_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()

# Install to package directory
install(TARGETS crm_cpp LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/..)
