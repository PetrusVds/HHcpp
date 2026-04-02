cmake_minimum_required(VERSION 3.15)
project(Hodgkin-Huxley)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

# --- Python Build Logic ---
if(SKBUILD) # Python package
    set(DEFAULT_BUILD_TESTS OFF)
    set(DEFAULT_BUILD_MAIN OFF)
else()
    set(DEFAULT_BUILD_TESTS ON)
    set(DEFAULT_BUILD_MAIN ON)
endif()

option(BUILD_TESTING "Build the C++ test suite" ${DEFAULT_BUILD_TESTS})
option(BUILD_MAIN "Build the standalone main executable" ${DEFAULT_BUILD_MAIN})
option(ENABLE_COVERAGE "Enable coverage reporting" OFF)

# --- Coverage Logic ---
if(ENABLE_COVERAGE)
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(--coverage)
        add_link_options(--coverage)
    endif()
endif()

# --- Core Library & Optimizations ---
file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
list(FILTER CORE_SOURCES EXCLUDE REGEX "main\\.cpp$")

add_library(project_core STATIC ${CORE_SOURCES})
set_target_properties(project_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(project_core PUBLIC include)

if(MSVC)
    target_compile_options(project_core INTERFACE $<$<CONFIG:Release,RelWithDebInfo>:/O2 /Ot>)
else()
    target_compile_options(project_core INTERFACE $<$<CONFIG:Release,RelWithDebInfo>:-O3>)
    target_compile_options(project_core INTERFACE $<$<CONFIG:Release,RelWithDebInfo>:-march=native>)
endif()

# --- Executables and Modules ---
if(BUILD_MAIN)
    add_executable(main src/main.cpp)
    target_link_libraries(main PRIVATE project_core)
endif()

find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    include(FetchContent)
    FetchContent_Declare(pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11 GIT_TAG v2.12.0)
    FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(HHcpp python/bindings.cpp)
target_link_libraries(HHcpp PRIVATE project_core)
install(TARGETS HHcpp DESTINATION .)

# --- Tests ---
if(BUILD_TESTING)
    enable_testing()
    find_package(Boost REQUIRED)
    
    file(GLOB TEST_FILES "tests/test-*.cpp")
    foreach(test_source ${TEST_FILES})
        get_filename_component(test_name ${test_source} NAME_WE)
        add_executable(${test_name} ${test_source})
        target_link_libraries(${test_name} PRIVATE project_core Boost::boost)
        add_test(NAME ${test_name} COMMAND ${test_name})    
    endforeach()
endif()