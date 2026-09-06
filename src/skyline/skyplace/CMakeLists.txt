add_subdirectory(src/gpu_louvain)

find_package(Eigen3 REQUIRED)
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core Widgets)
find_package(OpenMP REQUIRED)
find_package(CUDAToolkit)
enable_language(CUDA)

add_library(skyplace)

# CPU part
target_sources(skyplace
	PRIVATE
  src/SkyPlace.cpp
  src/SkyPlaceDB.cpp
  src/Util.cpp
  src/painter/Painter.cpp

  src/initial_place/InitialPlacer.cpp
  src/initial_place/FMPartitioner.cpp
  src/initial_place/SDPSolverCPU.cpp

  src/object/GPRect.cpp
  src/object/GPCell.cpp
  src/object/GPNet.cpp
  src/object/GPPin.cpp
  src/object/GPBin.cpp
)

# GPU part
target_sources(skyplace
  PRIVATE
  src/nonlinear_solver/SolverBase.cu
  src/nonlinear_solver/NesterovSolver.cu
  src/nonlinear_solver/AdamSolver.cu
  src/wirelength_gradient/WireLengthGradient.cu
  src/density_gradient/PoissonSolver.cu
  src/density_gradient/CudaDct.cu
  src/density_gradient/DensityGradient.cu
  src/initial_place/SDPSolverGPU.cu
  src/problem_instance/TargetFunction.cu
)

target_include_directories(skyplace
  PUBLIC
  include 
  PRIVATE
  src
  /home/jkim/mosek/10.1/tools/platform/linux64x86/h
)

# CPU stuffs
target_link_libraries(skyplace
  PRIVATE
  db
  cuda_linalg
  gpu_louvain 
	Eigen3::Eigen
  Qt${QT_VERSION_MAJOR}::Core
  Qt${QT_VERSION_MAJOR}::Widgets
  OpenMP::OpenMP_CXX
  /home/jkim/mosek/10.1/tools/platform/linux64x86/bin/libmosek64.so.10.1
  /home/jkim/mosek/10.1/tools/platform/linux64x86/bin/libfusion64.so.10.1
)

# GPU stuffs
target_link_libraries(skyplace
  PRIVATE
  CUDA::cufft
  cuda_linalg
)

set_target_properties(skyplace PROPERTIES CUDA_ARCHITECTURES "75")
