cmake_minimum_required(VERSION 3.17)
project(CRM_Dynamics VERSION 1.0.0 LANGUAGES CXX)

set(PROJECT_DESCRIPTION "Cosserat Rod Model for MRI-Actuated Robotic Catheter with ML/AI Integration")
set(PROJECT_HOMEPAGE_URL "https://github.com/yourusername/CRM_Dynamics")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(BUILD_PYTHON_BINDINGS "Build Python bindings with pybind11" ON)
option(BUILD_TESTING "Build tests" ON)
option(BUILD_EXAMPLES "Build example executables" ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option(ENABLE_NATIVE_ARCH "Enable native architecture optimizations" ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(ENABLE_NATIVE_ARCH AND CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found - enabling parallelization")
    else()
        message(WARNING "OpenMP not found - disabling parallelization")
        set(ENABLE_OPENMP OFF)
    endif()
endif()

if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter Development QUIET)
    find_package(pybind11 CONFIG QUIET)
    if(pybind11_FOUND AND Python3_FOUND)
        message(STATUS "Python ${Python3_VERSION} and pybind11 found - building Python bindings")
    else()
        message(WARNING "Python3 or pybind11 not found - skipping Python bindings")
        set(BUILD_PYTHON_BINDINGS OFF)
    endif()
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(NUM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/numerical")
set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/main")
set(PYTHON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python")

set(CRM_HEADERS
    ${SRC_DIR}/CRM.hpp
    ${SRC_DIR}/CRM_MatrixOperations.hpp
    ${SRC_DIR}/CRM_StateVector_Definitions.hpp
    ${SRC_DIR}/CRM_BVPIVP_APIDeclarations.hpp
    ${SRC_DIR}/CRM_IVP_NumericalIntegrationTemplates.hpp
    ${SRC_DIR}/CRM_IVPJacobian.hpp
    ${SRC_DIR}/CRM_IVPJacobian_InternalAPI.hpp
    ${SRC_DIR}/CRM_FK_InternalAPI.hpp
    ${SRC_DIR}/CoilDynamics_Defs.cpp
    ${SRC_DIR}/CRMDYN.hpp
    ${SRC_DIR}/CRMDYN_Numerical_Integration.hpp
)

set(CRM_SOURCES
    ${SRC_DIR}/CRM_IVPSolver.cpp
    ${SRC_DIR}/CRM_BVPSolver.cpp
    ${SRC_DIR}/CRM_IVPJacobian.cpp
    ${SRC_DIR}/CRM_ForwardKinematics.cpp
    ${SRC_DIR}/CRM_CatheterClass.cpp
    ${SRC_DIR}/CRM_SupportFunctions.cpp
)

set(NUMERICAL_HEADERS
    ${NUM_DIR}/minpack.hpp
    ${NUM_DIR}/minpack_TemplateDefinitions.hpp
    ${NUM_DIR}/minpack_wGivenJac_TemplateDefinitions.hpp
    ${NUM_DIR}/minpack_DYN.hpp
)

set(NUMERICAL_SOURCES
    ${NUM_DIR}/minpack.cpp
    ${NUM_DIR}/minpack_DYN_Defs.cpp
)

add_library(CRMCPPLib STATIC
    ${CRM_HEADERS}
    ${CRM_SOURCES}
    ${NUMERICAL_HEADERS}
    ${NUMERICAL_SOURCES}
)

target_include_directories(CRMCPPLib PUBLIC
    $<BUILD_INTERFACE:${SRC_DIR}>
    $<BUILD_INTERFACE:${NUM_DIR}>
    $<INSTALL_INTERFACE:include/CRM>
)

target_link_libraries(CRMCPPLib PUBLIC Eigen3::Eigen)

if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(CRMCPPLib PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(CRMCPPLib PUBLIC CRM_USE_OPENMP)
endif()

target_compile_options(CRMCPPLib PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g -O0>
)

set_target_properties(CRMCPPLib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${CRM_HEADERS}"
)

if(BUILD_EXAMPLES)
    add_executable(CRMTest
        ${MAIN_DIR}/CRMTest.cpp
        ${MAIN_DIR}/CRMTest.h
        ${MAIN_DIR}/CRM_KinematicsTestFunctions.hpp
        ${MAIN_DIR}/CRM_KinematicsTestFunctions.cpp
    )
    target_include_directories(CRMTest PRIVATE ${SRC_DIR} ${NUM_DIR})
    target_link_libraries(CRMTest PRIVATE CRMCPPLib)
    target_compile_options(CRMTest PRIVATE -Wall -Wextra)

    add_executable(CRMDYNTest
        ${MAIN_DIR}/CRMDYN_test.cpp
        ${MAIN_DIR}/CRMTest.h
        ${MAIN_DIR}/CRM_KinematicsTestFunctions.hpp
        ${MAIN_DIR}/CRM_KinematicsTestFunctions.cpp
    )
    target_include_directories(CRMDYNTest PRIVATE ${SRC_DIR} ${NUM_DIR})
    target_link_libraries(CRMDYNTest PRIVATE CRMCPPLib)
    target_compile_options(CRMDYNTest PRIVATE -Wall -Wextra)
endif()

if(BUILD_PYTHON_BINDINGS AND pybind11_FOUND AND Python3_FOUND)
    if(EXISTS "${PYTHON_DIR}/crm_bindings_complete.cpp")
        pybind11_add_module(crm_cpp MODULE
            ${PYTHON_DIR}/crm_bindings_complete.cpp
        )
        target_include_directories(crm_cpp PRIVATE ${SRC_DIR} ${NUM_DIR})
        target_link_libraries(crm_cpp PRIVATE CRMCPPLib Eigen3::Eigen)
        
        if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
            target_link_libraries(crm_cpp PRIVATE OpenMP::OpenMP_CXX)
        endif()
        
        target_compile_options(crm_cpp PRIVATE
            -Wall -Wextra
            $<$<CONFIG:Release>:-O3>
        )
        
        install(TARGETS crm_cpp
            LIBRARY DESTINATION ${Python3_SITELIB}
            COMPONENT python
        )
    else()
        message(WARNING "Python bindings source not found: ${PYTHON_DIR}/crm_bindings_complete.cpp")
    endif()
endif()

if(BUILD_TESTING)
    enable_testing()
    
    if(TARGET CRMTest)
        add_test(NAME CRMKinematicsTest COMMAND CRMTest)
        set_tests_properties(CRMKinematicsTest PROPERTIES LABELS "kinematics" TIMEOUT 300)
    endif()
    
    if(TARGET CRMDYNTest)
        add_test(NAME CRMDynamicsTest COMMAND CRMDYNTest)
        set_tests_properties(CRMDynamicsTest PROPERTIES LABELS "dynamics" TIMEOUT 600)
    endif()
    
    if(BUILD_PYTHON_BINDINGS AND TARGET crm_cpp)
        add_test(NAME PythonModuleImportTest
            COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.path.insert(0, '${CMAKE_CURRENT_BINARY_DIR}'); import crm_cpp"
        )
        set_tests_properties(PythonModuleImportTest PROPERTIES LABELS "python" TIMEOUT 30)
    endif()
endif()

include(GNUInstallDirs)

install(TARGETS CRMCPPLib
    EXPORT CRMDynamicsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/CRM
    COMPONENT library
)

install(DIRECTORY ${SRC_DIR}/ ${NUM_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/CRM
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
    COMPONENT headers
)

if(BUILD_EXAMPLES)
    install(TARGETS CRMTest CRMDYNTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT executables
    )
endif()

if(BUILD_PYTHON_BINDINGS)
    install(DIRECTORY ${PYTHON_DIR}/
        DESTINATION ${Python3_SITELIB}/crm_catheter
        FILES_MATCHING 
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE
            PATTERN "*.cpp" EXCLUDE
        COMPONENT python
    )
endif()

install(EXPORT CRMDynamicsTargets
    FILE CRMDynamicsTargets.cmake
    NAMESPACE CRMDynamics::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CRMDynamics
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CRMDynamicsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CRMDynamicsConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CRMDynamics
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/CRMDynamicsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CRMDynamicsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/CRMDynamicsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CRMDynamics
)

set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")

include(CPack)
