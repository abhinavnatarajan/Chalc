cmake_minimum_required(VERSION 3.20...3.26)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(DEFINED ENV{VCPKG_ROOT_DIR})
    message(STATUS "Using toolchain $ENV{VCPKG_ROOT_DIR}/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_ROOT_DIR}/scripts/buildsystems/vcpkg.cmake)
endif()

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})

find_package(CGAL REQUIRED)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE src/main.cpp)
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_core PRIVATE CGAL::CGAL)

install(TARGETS _core DESTINATION Chalc)