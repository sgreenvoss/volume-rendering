cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(HighlightPickedActor)

find_package(Kokkos 5.1 REQUIRED CONFIG)
find_package(VTK COMPONENTS
  RenderingCore
  RenderingOpenGL2
  IOLegacy
  InteractionStyle
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "VolumeRendering: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(VolumeRendering MACOSX_BUNDLE
  VolumeRendering.cxx
  astro.cxx
)

target_link_libraries(VolumeRendering PRIVATE
  ${VTK_LIBRARIES}
  Kokkos::kokkos
)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS VolumeRendering
  MODULES ${VTK_LIBRARIES}
)