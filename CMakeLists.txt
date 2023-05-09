cmake_minimum_required(VERSION 3.23)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

if (DEFINED ENV{CMAKE_USE_VCPKG})
    set(CMAKE_USE_VCPKG $ENV{CMAKE_USE_VCPKG})
else()
    set(CMAKE_USE_VCPKG ON)
endif()

if ($ENV{CI})
    message(STATUS "CI build detected.")
    set(ENV{VCPKG_INSTALLATION_ROOT} ${CMAKE_CURRENT_LIST_DIR}/vcpkg)
    message(STATUS "Setting VCPKG_INSTALLATION_ROOT to ${CMAKE_CURRENT_LIST_DIR}/vcpkg")
endif()

# If Vcpkg is available and not explicitly disabled then use it
if (DEFINED ENV{VCPKG_INSTALLATION_ROOT} AND CMAKE_USE_VCPKG)
    if (EXISTS $ENV{VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
        message(STATUS "Using Vcpkg toolchain: ${CMAKE_TOOLCHAIN_FILE}")
        set(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
    else ()
        message(FATAL_ERROR "Did not find vcpkg at specified location: $ENV{VCPKG_INSTALLATION_ROOT}")
    endif()
# for local builds where vcpkg is not available
elseif ((NOT DEFINED ENV{VCPKG_INSTALLATION_ROOT}) AND CMAKE_USE_VCPKG)
    message(STATUS "Fetching vcpkg. This could take some time...")
    include(FetchContent)
        FetchContent_Declare(vcpkg
            GIT_REPOSITORY https://github.com/microsoft/vcpkg/
            GIT_TAG 2022.09.27
        )
    FetchContent_MakeAvailable(vcpkg)
    message(STATUS "Fetching vcpkg -- done")
    message(STATUS "Using Vcpkg toolchain: ${CMAKE_TOOLCHAIN_FILE}")
    set(CMAKE_TOOLCHAIN_FILE "${vcpkg_SOURCE_DIR}/scripts/buildsystems/vcpkg.cmake")
else ()
    message(WARNING "Proceeding without Vcpkg.")
endif()

if (DEFINED SKBUILD)
    message(STATUS "Detected scikit-build environment.")
    message(STATUS "Project name: ${SKBUILD_PROJECT_NAME}")
    message(STATUS "Project version: ${SKBUILD_PROJECT_VERSION}")
endif()
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)
set(SUBMODULE_NAME core)

# Find dependencies
find_package(Eigen3 CONFIG REQUIRED)
if (Eigen3_FOUND)
    message(STATUS "Found Eigen3 headers in: ${EIGEN3_INCLUDE_DIR}")
endif()
find_package(CGAL CONFIG REQUIRED)
if (DEFINED CGAL_MODULES_DIR AND Eigen3_FOUND)
    include(${CGAL_MODULES_DIR}/CGAL_Eigen3_support.cmake)
    message(STATUS "CGAL has Eigen3 support.")
endif()

# ----------------------------- Build python api ------------------------------

# Option for building the project Python API, `setup.py` defines this variable
if(DEFINED SKBUILD)
    message(STATUS "Build extension module...")

    message(STATUS "Python executable: ${PYTHON_EXECUTABLE}")
    # call pybind11-config to obtain the root of the cmake package
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR_RAW)
    string(STRIP ${pybind11_DIR_RAW} pybind11_DIR)

    find_package(pybind11 CONFIG REQUIRED)

    # Create python extension module ###
    pybind11_add_module(${SUBMODULE_NAME} MODULE src/cpp/pybind11_bindings.cpp src/cpp/filtration.cpp src/cpp/chromatic.cpp)
    target_link_libraries(${SUBMODULE_NAME} PRIVATE Eigen3::Eigen CGAL::CGAL CGAL::Eigen3_support)
    target_compile_definitions(${SUBMODULE_NAME} PRIVATE VERSION_INFO=${PROJECT_VERSION})

    if (CMAKE_USE_VCPKG) 
        # List `vcpkg` dependencies from VCPKG_INSTALLATION_ROOT so we can install them
        # together with the extension module.
        file(GLOB RUNTIME_VCPKG_DEPENDENCIES
            RELATIVE "$ENV{VCPKG_INSTALLATION_ROOT}/${VCPKG_TARGET_TRIPLET}/bin/"
            "$ENV{VCPKG_INSTALLATION_ROOT}/${VCPKG_TARGET_TRIPLET}/bin/*")

        # Install the `vcpkg` dependencies that are used in our ${PROJECT_NAME}
        # library in the same DESTINATION as the python API. Note that it
        # only supports collecting the runtime dependencies for Windows, Linux and
        # macOS platforms.
        install(TARGETS ${SUBMODULE_NAME}
            RUNTIME_DEPENDENCIES
            PRE_INCLUDE_REGEXES ${RUNTIME_VCPKG_DEPENDENCIES}
            PRE_EXCLUDE_REGEXES ".*"
            DESTINATION .)
    else()
        # we assume that runtime dependencies are system libraries or otherwise available in LD_LIBRARY_PATH
        install(TARGETS ${SUBMODULE_NAME} DESTINATION .)
    endif()
endif()