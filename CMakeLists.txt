cmake_minimum_required(VERSION 3.23)
cmake_policy(SET CMP0054 NEW)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT DEFINED ENV{NO_USE_VCPKG})
    if (DEFINED ENV{VCPKG_INSTALLATION_ROOT})
        set(VCPKG_INSTALLATION_ROOT $ENV{VCPKG_INSTALLATION_ROOT})
    elseif ($ENV{CI})
        message(STATUS "CI build detected.")
        set(VCPKG_INSTALLATION_ROOT ${CMAKE_CURRENT_LIST_DIR}/vcpkg)
        message(STATUS "Setting VCPKG_INSTALLATION_ROOT to ${CMAKE_CURRENT_LIST_DIR}/vcpkg")
    else()
        message(STATUS "Downloading vcpkg from https://github.com/microsoft/vcpkg/")
        message(STATUS "Commit SHA512: ${VCPKG_COMMIT_ID}")
        include(FetchContent)
            FetchContent_Declare(vcpkg
                GIT_REPOSITORY https://github.com/microsoft/vcpkg/
                GIT_TAG VCPKG_COMMIT_ID
            )
        FetchContent_MakeAvailable(vcpkg)
        message(STATUS "Downloading vcpkg -- done")
        set(VCPKG_INSTALLATION_ROOT ${vcpkg_SOURCE_DIR})
        message(STATUS "Bootstrapping vcpkg")
        execute_process(COMMAND ${VCPKG_INSTALLATION_ROOT}/bootstrap-vcpkg.bat)
    endif()
    if (EXISTS ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
        set(X_VCPKG_APPLOCAL_DEPS_INSTALL ON CACHE BOOL "Automatically copy dependencies into the install target directory for executables.")
        message(STATUS "Using Vcpkg toolchain: ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake")
        set(CMAKE_TOOLCHAIN_FILE ${VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake)
    else ()
        message(FATAL_ERROR "Did not find vcpkg at specified location: ${VCPKG_INSTALLATION_ROOT}")
    endif()
else()
    message(STATUS "Proceeding without vcpkg.")
endif()

message(STATUS "Project name: ${SKBUILD_PROJECT_NAME}")
message(STATUS "Project version: ${SKBUILD_PROJECT_VERSION}")
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)

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
pybind11_add_module(
    filtration 
    MODULE 
        src/cpp/filtration_py.cxx 
        src/cpp/filtration.cxx)

pybind11_add_module(
    chromatic 
    MODULE 
        src/cpp/chromatic_py.cxx 
        src/cpp/chromatic.cxx 
        src/cpp/filtration.cxx)
target_include_directories(chromatic PRIVATE src/cpp/ConstrainedMiniball/mpreal)
target_link_libraries(chromatic PRIVATE Eigen3::Eigen CGAL::CGAL CGAL::Eigen3_support)

pybind11_add_module(version MODULE src/cpp/version.cxx)
target_compile_definitions(version PRIVATE VERSION_INFO=${PROJECT_VERSION})

if (CMAKE_USE_VCPKG AND WIN32)
    # file(GLOB VCPKG_SHARED_LIBRARIES
    #     RELATIVE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/"
    #     "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/*")
    # the following sets the variable VCPKG_SHARED_LIBRARIES
    # install(CODE "
    #     message(STATUS \"Searching for runtime dependencies in: \")
    #     message(STATUS \"${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin\")

    #     if (NOT EXISTS \"${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin\")
    #         message(FATAL_ERROR \"Error when searching for runtime dependencies.\")
    #     endif()
    #     get_directory_property(targets PROPERTY BUILDSYSTEM_TARGETS)
    #     set(resolved_deps \"\")
    #     foreach(tgt \$\{targets\})
    #         file(GET_RUNTIME_DEPENDENCIES
    #             RESOLVED_DEPENDENCIES_VAR temp
    #             MODULES $<TARGET_FILE:\$\{tgt\}>
    #             PRE_INCLUDE_REGEXES ${VCPKG_SHARED_LIBRARIES}
    #             PRE_EXCLUDE_REGEXES \".*\"
    #             DIRECTORIES \"${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin\")
    #         LIST(APPEND resolved_deps \$\{temp\})
    #     endforeach()
        
    #     message(STATUS \"The following runtime dependencies will be copied into the distribution:\")
    #     message(STATUS \$\{resolved_deps\})
    #     ")
    # Install the `vcpkg` dependencies that are used in our ${PROJECT_NAME}
    # library in the same DESTINATION as the python API. Note that it
    # only supports collecting the runtime dependencies for Windows, Linux and
    # macOS platforms.
    # install(TARGETS filtration chromatic version
    #     RUNTIME_DEPENDENCIES
    #     PRE_INCLUDE_REGEXES ${VCPKG_SHARED_LIBRARIES}
    #     PRE_EXCLUDE_REGEXES ".*"
    #     DIRECTORIES "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin"
    #     DESTINATION .)
#     install(TARGETS filtration chromatic version
#         RUNTIME DESTINATION .
#         LIBRARY DESTINATION .)
# else()
#     # we assume that runtime dependencies are system libraries or otherwise available in LD_LIBRARY_PATH
#     install(TARGETS ${PROJECT_MODULE_NAME} DESTINATION .)
endif()
install(TARGETS chromatic filtration version DESTINATION .)