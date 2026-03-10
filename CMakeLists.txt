cmake_minimum_required(VERSION 3.16)
project(CTSurface)

set(CMAKE_CXX_STANDARD 17)

find_package(VTK REQUIRED COMPONENTS
    CommonCore CommonDataModel IOImage FiltersCore FiltersSources
    RenderingCore RenderingOpenGL2 InteractionStyle IOXML IOGeometry
    ImagingCore ImagingStencil ImagingGeneral)

find_package(OpenCV REQUIRED)
find_package(OpenMP)

add_executable(CTSurface
    main.cpp
    VolumeData.h
    SurfaceDetermination.h
    MarchingCubes.h
    VoidAnalysis.h
    mc_tables.h)


target_link_libraries(CTSurface
    ${VTK_LIBRARIES}
    ${OpenCV_LIBS})

if(OpenMP_CXX_FOUND)
    target_link_libraries(CTSurface OpenMP::OpenMP_CXX)
endif()

vtk_module_autoinit(
    TARGETS CTSurface
    MODULES ${VTK_LIBRARIES})
