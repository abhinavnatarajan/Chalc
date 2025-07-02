cmake_minimum_required(VERSION 4.0.0)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
cmake_policy(SET CMP0167 NEW) # Let Boost find its own BoostConfig.cmake

if(NOT DEFINED ENV{NO_USE_VCPKG})
    if(DEFINED ENV{CI})
        message(STATUS "CI build detected.")
        set(VCPKG_INSTALLATION_ROOT ${CMAKE_CURRENT_LIST_DIR}/vcpkg)
        message(
            STATUS
            "Setting VCPKG_INSTALLATION_ROOT to ${CMAKE_CURRENT_LIST_DIR}/vcpkg"
        )
    elseif(DEFINED ENV{VCPKG_ROOT})
        set(VCPKG_INSTALLATION_ROOT $ENV{VCPKG_ROOT})
        message(STATUS "VCPKG_INSTALLATION_ROOT is ${VCPKG_INSTALLATION_ROOT}")
    else()
        # get the vcpkg commit id from vcpkg.json
        file(READ ${CMAKE_CURRENT_LIST_DIR}/vcpkg.json VCPKG_JSON_STRING)
        string(JSON VCPKG_COMMIT_ID GET ${VCPKG_JSON_STRING} builtin-baseline)
        message(
            STATUS
            "Downloading vcpkg from https://github.com/microsoft/vcpkg/"
        )
        message(STATUS "Commit SHA512: ${VCPKG_COMMIT_ID}")
        include(FetchContent)
        FetchContent_Declare(
            vcpkg
            GIT_REPOSITORY https://github.com/microsoft/vcpkg/
            GIT_TAG ${VCPKG_COMMIT_ID}
        )
        FetchContent_MakeAvailable(vcpkg)
        message(STATUS "Downloading vcpkg -- done")
        set(VCPKG_INSTALLATION_ROOT ${vcpkg_SOURCE_DIR})
        message(STATUS "Bootstrapping vcpkg")
        if(WIN32)
            execute_process(
                COMMAND ${VCPKG_INSTALLATION_ROOT}/bootstrap-vcpkg.bat
            )
        elseif(UNIX)
            execute_process(
                COMMAND ${VCPKG_INSTALLATION_ROOT}/bootstrap-vcpkg.sh
            )
        endif()
    endif()
    if(EXISTS ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
        message(
            STATUS
            "Using Vcpkg toolchain: ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake"
        )
        set(CMAKE_TOOLCHAIN_FILE
            "${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake"
            CACHE STRING
            ""
        )
        set(USING_VCPKG TRUE)
    else()
        message(
            FATAL_ERROR
            "Did not find vcpkg at specified location: ${VCPKG_INSTALLATION_ROOT}"
        )
    endif()
else()
    message(STATUS "Proceeding without vcpkg.")
endif()

message(STATUS "Project name: ${SKBUILD_PROJECT_NAME}")
message(STATUS "Project version: ${SKBUILD_PROJECT_VERSION}")
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# Find dependencies
find_package(Eigen3 CONFIG REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Found Eigen3 headers in: ${EIGEN3_INCLUDE_DIR}")
endif()

find_package(TBB REQUIRED)

find_package(CGAL CONFIG REQUIRED)
if(DEFINED CGAL_MODULES_DIR AND Eigen3_FOUND)
    include(${CGAL_MODULES_DIR}/CGAL_Eigen3_support.cmake)
    message(STATUS "CGAL has Eigen3 support.")
endif()
if(DEFINED CGAL_MODULES_DIR AND TBB_FOUND)
    include(${CGAL_MODULES_DIR}/CGAL_TBB_support.cmake)
    message(STATUS "CGAL has TBB support.")
endif()

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# ----------------------------- Build python api ------------------------------

add_library(_filtration STATIC src/chalc/filtration/filtration.cxx)
target_include_directories(_filtration PUBLIC src)
add_library(_chromatic STATIC src/chalc/chromatic/chromatic.cxx)
target_link_libraries(
    _chromatic
	PUBLIC _filtration Eigen3::Eigen TBB::tbb CGAL::CGAL CGAL::Eigen3_support CGAL::TBB_support
)
set_target_properties(
    _filtration
    _chromatic
    PROPERTIES POSITION_INDEPENDENT_CODE ON
)

# Create python extension module
pybind11_add_module(filtration MODULE src/chalc/filtration_py.cxx)
target_link_libraries(filtration PRIVATE _filtration)
pybind11_add_module(chromatic MODULE src/chalc/chromatic_py.cxx)
target_link_libraries(chromatic PRIVATE _chromatic)

install(TARGETS filtration RUNTIME_DEPENDENCY_SET vcpkg_libs DESTINATION .)
install(TARGETS chromatic RUNTIME_DEPENDENCY_SET vcpkg_libs DESTINATION .)
if(USING_VCPKG AND WIN32)
    file(
        GLOB VCPKG_SHARED_LIBRARIES
        RELATIVE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/"
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/*"
    )
    install(
        RUNTIME_DEPENDENCY_SET
        vcpkg_libs
        DESTINATION
        lib
        PRE_INCLUDE_REGEXES ${VCPKG_SHARED_LIBRARIES}
        PRE_EXCLUDE_REGEXES ".*"
        DIRECTORIES "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/"
    )
endif()
