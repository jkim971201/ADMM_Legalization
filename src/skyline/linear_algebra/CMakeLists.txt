find_package(CUDAToolkit REQUIRED)
find_package(CCCL REQUIRED)
find_package(Eigen3 REQUIRED CONFIG)

enable_language(CUDA)

set(CUDA_LINALG_SRC
  src/CudaObject.cu
  src/CudaVector.cu
  src/CudaScalar.cu
  src/CudaVectorAlgebra.cu
  src/CudaMatrixAlgebra.cu
  src/CudaSparseMatrixAlgebra.cu
)

set_source_files_properties(${CUDA_LINALG_SRC} PROPERTIES LANGUAGE CUDA)

add_library(cuda_linalg OBJECT ${CUDA_LINALG_SRC})

#set_target_properties(cuda_linalg PROPERTIES CUDA_ARCHITECTURES 86)

target_link_libraries(cuda_linalg
  PRIVATE
  Eigen3::Eigen

  PUBLIC
  CUDA::cudart
  CUDA::cublas
  CUDA::cusparse
  CUDA::cusolver
  CCCL::Thrust
)

target_include_directories(cuda_linalg
  PUBLIC
  include
  PRIVATE
  Eigen3::Eigen
)
