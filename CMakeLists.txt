cmake_minimum_required(VERSION 3.23)
cmake_policy(SET CMP0054 NEW)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

if (DEFINED ENV{CMAKE_USE_VCPKG})
    set(CMAKE_USE_VCPKG $ENV{CMAKE_USE_VCPKG})
else()
    set(CMAKE_USE_VCPKG ON)
endif()

if (DEFINED ENV{VCPKG_INSTALLATION_ROOT})
    set(VCPKG_INSTALLATION_ROOT $ENV{VCPKG_INSTALLATION_ROOT})
endif()

if ($ENV{CI})
    message(STATUS "CI build detected.")
    set(VCPKG_INSTALLATION_ROOT ${CMAKE_CURRENT_LIST_DIR}/vcpkg)
    message(STATUS "Setting VCPKG_INSTALLATION_ROOT to ${CMAKE_CURRENT_LIST_DIR}/vcpkg")
endif()

# If Vcpkg is available and not explicitly disabled then use it
if (DEFINED VCPKG_INSTALLATION_ROOT AND CMAKE_USE_VCPKG)
    if (EXISTS ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
        message(STATUS "Using Vcpkg toolchain: ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake")
        set(CMAKE_TOOLCHAIN_FILE ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
    else ()
        message(FATAL_ERROR "Did not find vcpkg at specified location: ${VCPKG_INSTALLATION_ROOT}")
    endif()
# for local builds where vcpkg is not available
elseif ((NOT DEFINED VCPKG_INSTALLATION_ROOT) AND CMAKE_USE_VCPKG)
    message(STATUS "Fetching vcpkg. This could take some time...")
    include(FetchContent)
        FetchContent_Declare(vcpkg
            GIT_REPOSITORY https://github.com/microsoft/vcpkg/
            GIT_TAG VCPKG_COMMIT_ID
        )
    FetchContent_MakeAvailable(vcpkg)
    message(STATUS "Fetching vcpkg -- done")
    message(STATUS "Using Vcpkg toolchain: ${vcpkg_SOURCE_DIR}/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "${vcpkg_SOURCE_DIR}/scripts/buildsystems/vcpkg.cmake")
    set(VCPKG_INSTALLATION_ROOT ${vcpkg_SOURCE_DIR})
else ()
    message(WARNING "Proceeding without Vcpkg.")
endif()

message(STATUS "Project name: ${SKBUILD_PROJECT_NAME}")
message(STATUS "Project version: ${SKBUILD_PROJECT_VERSION}")
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)
set(PROJECT_MODULE_NAME core)

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

message(STATUS "Building extension module...")

message(STATUS "Python executable: ${PYTHON_EXECUTABLE}")
# call pybind11-config to obtain the root of the cmake package
execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_DIR_RAW)
string(STRIP ${pybind11_DIR_RAW} pybind11_DIR)

find_package(pybind11 CONFIG REQUIRED)

# Create python extension module ###
pybind11_add_module(${PROJECT_MODULE_NAME} MODULE src/cpp/pybind11_bindings.cpp src/cpp/filtration.cpp src/cpp/chromatic.cpp)
target_link_libraries(${PROJECT_MODULE_NAME} PRIVATE Eigen3::Eigen CGAL::CGAL CGAL::Eigen3_support)
target_compile_definitions(${PROJECT_MODULE_NAME} PRIVATE VERSION_INFO=${PROJECT_VERSION})

if (CMAKE_USE_VCPKG AND WIN32)
    file(GLOB VCPKG_SHARED_LIBRARIES
        RELATIVE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/"
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/*")
    # the following sets the variable VCPKG_SHARED_LIBRARIES
    install(CODE "
        message(STATUS \"Searching for runtime dependencies in: \")
        message(STATUS \"${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin\")

        if (NOT EXISTS \"${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin\")
            message(FATAL_ERROR \"Error when searching for runtime dependencies.\")
        endif()

        file(GET_RUNTIME_DEPENDENCIES
            RESOLVED_DEPENDENCIES_VAR resolved_deps
            MODULES $<TARGET_FILE:${PROJECT_MODULE_NAME}>
            PRE_INCLUDE_REGEXES ${VCPKG_SHARED_LIBRARIES}
            PRE_EXCLUDE_REGEXES \".*\"
            DIRECTORIES \"${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin\")
        
        message(STATUS \"The following runtime dependencies will be copied into the distribution:\")
        message(STATUS \$\{resolved_deps\})
        ")
    # Install the `vcpkg` dependencies that are used in our ${PROJECT_NAME}
    # library in the same DESTINATION as the python API. Note that it
    # only supports collecting the runtime dependencies for Windows, Linux and
    # macOS platforms.
    install(TARGETS ${PROJECT_MODULE_NAME}
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_SHARED_LIBRARIES}
        PRE_EXCLUDE_REGEXES ".*"
        DIRECTORIES "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin"
        DESTINATION .)
else()
    # we assume that runtime dependencies are system libraries or otherwise available in LD_LIBRARY_PATH
    install(TARGETS ${PROJECT_MODULE_NAME} DESTINATION .)
endif()