cmake_minimum_required(VERSION 3.17)
project(CRM_Dynamics)

# ============================================================
# C++ STANDARD & BUILD CONFIGURATION
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# ============================================================
# POSITION INDEPENDENT CODE (for shared libraries/Python)
# ============================================================
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")

# ============================================================
# FIND REQUIRED PACKAGES
# ============================================================

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Python3 COMPONENTS Interpreter Development QUIET)
find_package(pybind11 CONFIG QUIET)

# ============================================================
# PROJECT STRUCTURE SETUP
# ============================================================

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

# Source directory variables
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(NUM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/numerical")
set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/main")

# ============================================================
# CORE COSSERAT ROD MODEL SOURCES
# ============================================================

set(CRMSOURCES
    ${SRC_DIR}/CRM.hpp
    ${SRC_DIR}/CRM_MatrixOperations.hpp
    ${SRC_DIR}/CRM_StateVector_Definitions.hpp
    ${SRC_DIR}/CRM_BVPIVP_APIDeclarations.hpp
    ${SRC_DIR}/CRM_IVP_NumericalIntegrationTemplates.hpp
    ${SRC_DIR}/CRM_IVPSolver.cpp
    ${SRC_DIR}/CRM_BVPSolver.cpp
    ${SRC_DIR}/CRM_IVPJacobian.hpp
    ${SRC_DIR}/CRM_IVPJacobian_InternalAPI.hpp
    ${SRC_DIR}/CRM_IVPJacobian.cpp
    ${SRC_DIR}/CRM_FK_InternalAPI.hpp
    ${SRC_DIR}/CRM_ForwardKinematics.cpp
    ${SRC_DIR}/CRM_CatheterClass.cpp
    ${SRC_DIR}/CRM_SupportFunctions.cpp
    ${SRC_DIR}/CoilDynamics_Defs.cpp
    ${SRC_DIR}/CRMDYN.hpp
    ${SRC_DIR}/CRMDYN_Numerical_Integration.hpp
    ${NUM_DIR}/minpack.hpp
    ${NUM_DIR}/minpack_TemplateDefinitions.hpp
    ${NUM_DIR}/minpack_wGivenJac_TemplateDefinitions.hpp
    ${NUM_DIR}/minpack.cpp
    ${NUM_DIR}/minpack_DYN.hpp
    ${NUM_DIR}/minpack_DYN_Defs.cpp
)

# ============================================================
# BUILD: STATIC LIBRARY
# ============================================================

add_library(CRMCPPLib STATIC ${CRMSOURCES})

target_include_directories(CRMCPPLib PUBLIC
    "${SRC_DIR}"
    "${NUM_DIR}"
)

target_link_libraries(CRMCPPLib PUBLIC Eigen3::Eigen)

target_compile_options(CRMCPPLib PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3 -march=native>
    $<$<CONFIG:Debug>:-g -O0>
)

# ============================================================
# BUILD: EXECUTABLE - CRMTest (Kinematics)
# ============================================================

add_executable(CRMTest
    ${MAIN_DIR}/CRMTest.cpp
    ${MAIN_DIR}/CRMTest.h
    ${MAIN_DIR}/CRM_KinematicsTestFunctions.hpp
    ${MAIN_DIR}/CRM_KinematicsTestFunctions.cpp
)

target_include_directories(CRMTest PUBLIC
    "${SRC_DIR}"
    "${NUM_DIR}"
)

target_link_libraries(CRMTest PUBLIC
    CRMCPPLib
    Eigen3::Eigen
)

target_compile_options(CRMTest PRIVATE
    -Wall -Wextra
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g>
)

# ============================================================
# BUILD: EXECUTABLE - CRMDYNTest (Dynamics)
# ============================================================

add_executable(CRMDYNTest
    ${MAIN_DIR}/CRMDYN_test.cpp
    ${MAIN_DIR}/CRMTest.h
    ${MAIN_DIR}/CRM_KinematicsTestFunctions.hpp
    ${MAIN_DIR}/CRM_KinematicsTestFunctions.cpp
)

target_include_directories(CRMDYNTest PUBLIC
    "${SRC_DIR}"
    "${NUM_DIR}"
)

target_link_libraries(CRMDYNTest PUBLIC
    CRMCPPLib
    Eigen3::Eigen
)

target_compile_options(CRMDYNTest PRIVATE
    -Wall -Wextra
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g>
)

# ============================================================
# BUILD: PYTHON BINDINGS (pybind11)
# ============================================================

if(pybind11_FOUND AND Python3_FOUND)
    message(STATUS "Building Python bindings with pybind11")
    
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/python/crm_bindings.cpp")
        pybind11_add_module(crm_cpp
            python/crm_bindings.cpp
        )
        
        target_include_directories(crm_cpp PRIVATE
            "${SRC_DIR}"
            "${NUM_DIR}"
        )
        
        target_link_libraries(crm_cpp PRIVATE
            CRMCPPLib
            Eigen3::Eigen
        )
        
        target_compile_options(crm_cpp PRIVATE
            -Wall -Wextra
            $<$<CONFIG:Release>:-O3 -march=native>
        )
        
        install(TARGETS crm_cpp
            LIBRARY DESTINATION "${Python3_SITELIB}"
        )
    else()
        message(WARNING "python/crm_bindings.cpp not found; skipping Python bindings")
    endif()
else()
    message(WARNING "pybind11 or Python3 not found; skipping Python bindings")
endif()

# ============================================================
# INSTALLATION TARGETS
# ============================================================

install(TARGETS CRMTest CRMDYNTest
    RUNTIME DESTINATION bin
)

install(TARGETS CRMCPPLib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ${SRC_DIR}/ DESTINATION include/CRM
    FILES_MATCHING PATTERN "*.hpp"
)

# ============================================================
# TESTING & CPACK
# ============================================================

enable_testing()

add_test(NAME CRMTest COMMAND CRMTest)
add_test(NAME CRMDYNTest COMMAND CRMDYNTest)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION "1.0.0")
include(CPack)

# ============================================================
# BUILD SUMMARY
# ============================================================

message(STATUS "")
message(STATUS "========== Build Summary ==========")
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Position Independent Code: ON")
message(STATUS "Eigen3 Found: YES")
message(STATUS "Python3 Found: ${Python3_FOUND}")
message(STATUS "pybind11 Found: ${pybind11_FOUND}")
message(STATUS "====================================")
message(STATUS "")