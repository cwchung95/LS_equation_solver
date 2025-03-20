# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10...3.24)

# Project name
project(LS_Equation_Solver)

# Set the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(LIB_DIR ${CMAKE_SOURCE_DIR}/lib)
set(RUN_DIR ${CMAKE_SOURCE_DIR}/run)
set(OBJ_DIR ${CMAKE_SOURCE_DIR}/obj)

# Include directories
include_directories(${LIB_DIR})

# Retrive GSL flags using gsl-config
execute_process(COMMAND gsl-config --cflags OUTPUT_VARIABLE GSL_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND gsl-config --libs OUTPUT_VARIABLE GSL_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

# Add GSL flags to the project
add_compile_options(${GSL_CFLAGS})
link_libraries(${GSL_LIBS})

# Find ROOT package
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# Add warnings and errors
# add_compile_options(-Wall -Wextra -Wpedantic )

# Find all .cpp files in the src directory
file(GLOB SOURCES ${SRC_DIR}/*.cpp)

# Identify library source files based on matching headers in lib/
file(GLOB LIB_HEADERS ${LIB_HEADERS}/*.h ${LIB_DIR}/*.hpp)
set(LIB_SOURCES "")
foreach(header ${LIB_HEADERS})
  get_filename_component(header_base ${header} NAME_WE)
  list(APPEND LIB_SOURCES ${SRC_DIR}/${header_base}.cpp)
endforeach()

# Filter out library sources from the main sources
list(REMOVE_ITEM SOURCES ${LIB_SOURCES})

# Generate executables from non-library source files
# Modify the executable generation loop
foreach(source ${SOURCES})
    get_filename_component(exec_name ${source} NAME_WE)
    add_executable(${exec_name} ${source} ${LIB_SOURCES})
    set_target_properties(${exec_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${RUN_DIR})
    target_link_libraries(${exec_name} ${ROOT_LIBRARIES})
endforeach()

