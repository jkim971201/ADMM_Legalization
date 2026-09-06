find_package(CUDAToolkit REQUIRED)
find_package(Thrust)

enable_language(CUDA)

set(GPU_LOUVAIN_SRC
  utils.cu
	community_aggregation.cu
	modularity_optimisation.cu
	removeSparse.cu
	gpulouvain.cu
)

add_library(gpu_louvain OBJECT ${GPU_LOUVAIN_SRC})

set_source_files_properties(${GPU_LOUVAIN_SRC} PROPERTIES LANGUAGE CUDA)

set_target_properties(gpu_louvain	PROPERTIES CUDA_ARCHITECTURES "75")

target_link_libraries(gpu_louvain
	PUBLIC
	CUDA::cufft
	CUDA::cudart
)

target_include_directories(gpu_louvain
  PRIVATE
	${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
)

