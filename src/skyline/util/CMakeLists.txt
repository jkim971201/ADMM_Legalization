add_library(util)

target_sources(util
  PRIVATE
  src/Chrono.cpp
)

target_include_directories(util
  PUBLIC
  include
)
