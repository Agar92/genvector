cmake_minimum_required(VERSION 3.6)

project(genvector CXX C) # C is for the IDE

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ ISO Standard")
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h ${CMAKE_CURRENT_SOURCE_DIR}/inc/GenVector/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)

add_library(genvector ${sources} ${headers})
target_link_libraries(genvector ${OpenMP_CXX_LIBRARIES})
target_include_directories(genvector PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc/> $<INSTALL_INTERFACE:include>)
set_property(TARGET genvector PROPERTY CXX_STANDARD 17)


set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
install(TARGETS genvector DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/inc/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/)

enable_testing()
add_subdirectory(test)

find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_INPUT "${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/inc/Math/ ${CMAKE_CURRENT_SOURCE_DIR}/inc/Math/GenVector/")
    get_property(PROJECT_INCLUDE_DIRS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#     FIXME empty
    foreach(dir ${PROJECT_INCLUDE_DIRS})
      set(DOXYGEN_INCLUDE_DIRS "${DOXYGEN_INCLUDE_DIRS} ${dir}")
    endforeach()
    foreach(dir ${CMAKE_EXTRA_GENERATOR_CXX_SYSTEM_INCLUDE_DIRS})
      set(DOXYGEN_INCLUDE_DIRS "${DOXYGEN_INCLUDE_DIRS} ${dir}")
    endforeach()
    foreach(dir ${CMAKE_EXTRA_GENERATOR_C_SYSTEM_INCLUDE_DIRS})
      set(DOXYGEN_INCLUDE_DIRS "${DOXYGEN_INCLUDE_DIRS} ${dir}")
    endforeach()
    # set input and output files
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message(STATUS "Configured doxygen")

#     message(STATUS "INPUT for doxygen is ${DOXYGEN_INPUT}")
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

