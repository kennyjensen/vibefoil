cmake_minimum_required(VERSION 3.20)

project(xfoil_wasm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".js")
endif()

add_library(xfoil_core
    naca.cpp
    spline.cpp
    userio.cpp
    xbl.cpp
    xblcom.cpp
    xblsys.cpp
    xfoil.cpp
    xgdes.cpp
    xgeom.cpp
    xoper.cpp
    xpanel.cpp
    xsolve.cpp
    xutils.cpp
)

target_include_directories(xfoil_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(xfoil_core PUBLIC cxx_std_17)

add_executable(compare_xbl
    compare_xbl.cpp
)

target_link_libraries(compare_xbl PRIVATE xfoil_core)

if(EMSCRIPTEN)
    target_link_options(compare_xbl PRIVATE -sWASM=1 -sALLOW_MEMORY_GROWTH=1 -sNODERAWFS=1)
endif()
