add_subdirectory(util)
add_subdirectory(linear_algebra)
add_subdirectory(database)
add_subdirectory(skygui)
add_subdirectory(legalizer)
add_subdirectory(skyplace)

add_library(skyline)

target_sources(skyline
  PRIVATE
  SkyLine.cpp
)

target_include_directories(skyline
  PUBLIC
	include
)

target_link_libraries(skyline
  db
	gui
  skyplace
  legalizer
  util
)
