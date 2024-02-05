cmake_minimum_required(VERSION 3.17)
project(CRM_Dynamics)
find_package(Eigen3 REQUIRED NO_MODULE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/../
        )

set(CRMSOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_MatrixOperations.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_StateVector_Definitions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_BVPIVP_APIDeclarations.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVP_NumericalIntegrationTemplates.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVPSolver.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_BVPSolver.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVPJacobian.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVPJacobian_InternalAPI.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVPJacobian.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_FK_InternalAPI.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_ForwardKinematics.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_CatheterClass.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_SupportFunctions.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_TemplateDefinitions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_wGivenJac_TemplateDefinitions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_DYN.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_DYN_Defs.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/src/CoilDynamics_Defs.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/src/CRMDYN.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/src/CRMDYN_Numerical_Integration.hpp
		)

add_library(CRMCPPLib STATIC
	${CRMSOURCES}
	)

target_include_directories(CRMCPPLib PUBLIC
	"${PROJECT_SOURCE_DIR}/src"
	"${PROJECT_SOURCE_DIR}/numerical"
        )

add_executable(CRMDYNTest
	main/CRMDYN_test.cpp
	main/CRMTest.h
	main/CRM_KinematicsTestFunctions.hpp
	main/CRM_KinematicsTestFunctions.cpp
	${CRMSOURCES}
	)

add_executable(CRMTest
		main/CRMTest.cpp
		main/CRMTest.h
		main/CRM_KinematicsTestFunctions.hpp
		main/CRM_KinematicsTestFunctions.cpp
		${CRMSOURCES}
		)

target_link_libraries(CRMTest PUBLIC
	Eigen3::Eigen
	)

target_link_libraries(CRMDYNTest PUBLIC
		Eigen3::Eigen
		)

target_include_directories(CRMTest PUBLIC
	"${PROJECT_SOURCE_DIR}/src"
	"${PROJECT_SOURCE_DIR}/numerical"
        )

target_include_directories(CRMDYNTest PUBLIC
		"${PROJECT_SOURCE_DIR}/src"
		"${PROJECT_SOURCE_DIR}/numerical"
		)

