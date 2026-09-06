cmake_minimum_required(VERSION 3.26)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

project(SkyLine LANGUAGES CXX)

set(CMAKE_C_COMPILER "g++")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# For CUDA 20
set(CMAKE_CUDA_STANDARD 20) 
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 75)

# Home Directory
add_subdirectory(${PROJECT_SOURCE_DIR}/src/skyline)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tclCmd)

# Source Code
set(PROJECT_SRC 
	src/main.cpp
)

# Executable
add_executable(${PROJECT_NAME} ${PROJECT_SRC})

# Compile Definition
target_compile_definitions(${PROJECT_NAME}
  PUBLIC 
)

# Include Directory
target_include_directories(${PROJECT_NAME}
  PUBLIC
	${TCL_INCLUDE_PATH} 
)

# Link Library
target_link_libraries(${PROJECT_NAME} 
  PUBLIC
  ${TCL_LIBRARY}
	skyline
  tcl_cmd	
)
