cmake_minimum_required(VERSION 3.18)
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
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_BVPIVP_APIDeclarations.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_BVP_Definitions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVP_Definitions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVPJacobian_Declarations.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_IVPJacobian_Definitions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_FK_InternalAPIDeclarations.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/CRM_FK_Definitions.hpp
#        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_Declarations.hpp
#        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_Definitions.hpp
#        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_wGivenJac_Definitions.hpp
#        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_DYN.hpp
#        ${CMAKE_CURRENT_SOURCE_DIR}/numerical/minpack_DYN_Defs.hpp
	)

add_executable(CRMTest
	main/CRMTest.cpp
	main/CRMTest.h
	${CRMSOURCES}
	)

target_link_libraries(CRMTest PUBLIC
	Eigen3::Eigen
	)

target_include_directories(CRMTest PUBLIC
	"${PROJECT_SOURCE_DIR}/src"
	"${PROJECT_SOURCE_DIR}/src/numerical"
        )

