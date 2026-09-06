cmake_minimum_required(VERSION 3.10)

# Source Code
add_library(legalizer)

find_package(OpenMP REQUIRED)
find_package(CUDAToolkit)

enable_language(CUDA)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math") # use __expf
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${OpenMP_CXX_FLAGS}") # use OpenMP in .cu file

set(LG_SRC
  src/Legalizer.cpp
  src/DatabaseIO.cpp
  src/Util.cpp
  src/Grid.cpp
  src/Preplace.cpp
  src/objects/LGCell.cpp
  src/objects/LGCellTechInfo.cpp
  src/objects/LGRect.cpp

  # CUDA Codes
  src/database/CudaDatabase.cu
  src/admm_legalizer/CudaAdmmLegalizer.cu
  src/post_process/CudaPostProcess.cu
  src/post_process/CudaQpSolver.cu
)

target_sources(legalizer
  PRIVATE
  ${LG_SRC}
)

# Include Directory
target_include_directories(legalizer
  PUBLIC
  include 
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link Library
target_link_libraries(legalizer
  PRIVATE
  db
  util
  cuda_linalg
  OpenMP::OpenMP_CXX
)
