cmake_minimum_required(VERSION 3.19)
project(CRM_Dynamics)

find_package(PythonLibs REQUIRED)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)


set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/../
        )

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

#add_library(matrixOperation
#        matrixOperation/matrixOperation.cpp
#        )



#install ILQR as a lib
#add_library(ILQR ILQR/ILQR.cpp)

#add_executable(planning_main planning_main.cpp)

#target_link_libraries(ILQR ${PYTHON_LIBRARIES} CRM_Model::CRM Load_LAA_model::Load_LAA_model)

target_link_libraries(planning_main )

#target_include_directories(ILQR
#        INTERFACE
#        $<INSTALL_INTERFACE:include>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#        )
#
#install(TARGETS ILQR
#        EXPORT ILQR_export
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib)
#
#install(EXPORT ILQR_export
#        FILE
#        ILQRTargets.cmake
#        NAMESPACE
#        ILQR::
#        DESTINATION
#        lib/cmake/ILQR)
#
#install(FILES
#        ${CMAKE_CURRENT_SOURCE_DIR}/ILQR.h
#        ${CMAKE_CURRENT_SOURCE_DIR}/matplotlibcpp.h
#        DESTINATION
#        "include/ILQR")
#
#install(FILES
#        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ILQRConfig.cmake
#        DESTINATION "lib/cmake/ILQR")

