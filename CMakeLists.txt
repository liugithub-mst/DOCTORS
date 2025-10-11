cmake_minimum_required(VERSION 3.18) # CUDA language support started in 3.8, improved later

# Project: enalbe C++ adn CUDA
project(DOCTORS
    VERSION 2.0
    LANGUAGES CXX CUDA
)

# ----- Options ------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_WARNINGS "Enable extra compiler warnings" ON)


# ----- Standards Use C++14 or newer -----
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)

# Default to Release if build type not set (nice for VS code)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
endif()

# Enable all CUDA architectures (or pick one)
# Let users override via -DCMAKE_CUDA_ARCHITECHTURES=...
# 'native' requrie CMake 3.24+ and CUDA 11.8+; otherwiase set a sane list.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # set(CMAKE_CUDA_ARCHITECTUES native)   # auto-detect GPU on the build machine
    set(CMAKE_CUDA_ARCHITECTURES 86)  # For RTX 20xx (75 Turing) and 30xx/40xx (86 Ampere)
endif()

# ----- Warings (C++ only; CUDA uses host compiler flags automatically) -----
if(ENABLE_WARNINGS)
    if(MSVC)
        add_compile_options(/W4 /permissive-)
    else()
        add_compile_options(-w -fexceptions)
    endif()
endif()

# ----- Directories -----
# Assume a layout like:
#   include/   # public headers
#   src/
#     main.cpp
#     lib/foo.cpp
#     lib/kernels.cu
#
set(SRC_DIR     "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Sanity checks (help catch path typos)
if(NOT EXISTS "${SRC_DIR}")
    message(FATAL_ERROR "SRC_DIR does not exist: ${SRC_DIR}")
endif()
if(NOT EXISTS "${INCLUDE_DIR}")
    message(WARNING "INCLUDE_DIR not found: ${INCLUDE_DIR}")
endif()

# Grab all .cpp and .cu under src/ (re-run cmake on file changes)
file(GLOB_RECURSE APP_SOURCES CONFIGURE_DEPENDS
    "${SRC_DIR}/*.cpp"
    "${SRC_DIR}/*.cu"
)

# ----- Executable -----
add_executable(${PROJECT_NAME} ${APP_SOURCES})

# Public headers live in include/; include them like #include "foo.cpp"
target_include_directories(${PROJECT_NAME}
    PRIVATE
        "${INCLUDE_DIR}"
)

# Optional: Enable separable compilation (for multiple .cu files)
set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link extra libraries if needed (like cudart)
target_link_libraries(${PROJECT_NAME} PRIVATE cuda cudart)
 