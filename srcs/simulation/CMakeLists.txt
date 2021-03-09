find_path(EIGEN_INCLUDE_DIRS NAMES "Eigen")
if(NOT EIGEN_INCLUDE_DIRS)
    message(FATAL_ERROR "Eigen not found!")
else()
	message(STATUS "Found Eigen: " ${EIGEN_INCLUDE_DIRS})
endif()

find_package(Lua 5.1 REQUIRED)

if(USE_MKL_WITH_TBB)
    find_package(MKL REQUIRED)
    if(MKL_FOUND)
        message(STATUS "MKL Found: " ${MKL_INCLUDE_DIR})
        message(STATUS "MKL Found: " ${MKL_CORE_LIB_ST})
        message(STATUS "MKL Found: " ${MKL_INTERFACE_INTEL_32BIT_LIB_ST})
        message(STATUS "MKL Found: " ${MKL_INTERFACE_INTEL_64BIT_LIB_ST})
        message(STATUS "MKL Found: " ${MKL_SEQ_LIB_ST})
        message(STATUS "MKL Found: " ${MKL_OMP_LIB_ST})
        message(STATUS "MKL Found: " ${MKL_TBB_LIB_ST})
    endif()
endif()

find_path(SOL3_INCLUDE_DIRS NAMES "sol/sol.hpp")
if(NOT SOL3_INCLUDE_DIRS)
    message(FATAL_ERROR "sol not found!")
else()
	message(STATUS "Found sol: " ${SOL3_INCLUDE_DIRS})
endif()

message(STATUS ${PROJECT_SOURCE_DIR})


file(GLOB_RECURSE LIBSIMULATION_SRCS
     "${PROJECT_SOURCE_DIR}/srcs/simulation/*.h"
     "${PROJECT_SOURCE_DIR}/srcs/simulation/*.inl"
     "${PROJECT_SOURCE_DIR}/srcs/simulation/*.hpp"
     "${PROJECT_SOURCE_DIR}/srcs/simulation/*.cpp")

add_library(pfemSimulation SHARED ${LIBSIMULATION_SRCS})
target_include_directories(pfemSimulation SYSTEM
                           PRIVATE ${GMSH_INCLUDE_DIRS}
                           PUBLIC ${EIGEN_INCLUDE_DIRS}
                           PUBLIC ${LUA_INCLUDE_DIR}
                           PUBLIC ${SOL3_INCLUDE_DIRS})
target_link_libraries(pfemSimulation
                      PRIVATE ${GMSH_LIBRARIES}
                      PUBLIC OpenMP::OpenMP_CXX
                      PUBLIC ${LUA_LIBRARIES} pfemMesh)
if(USE_MKL_WITH_TBB)
    target_link_libraries(pfemSimulation
                          PRIVATE mkl::mkl_intel_32bit_omp_dyn)
    target_compile_definitions(pfemSimulation PRIVATE EIGEN_USE_MKL_ALL)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    target_compile_options(pfemSimulation PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    target_compile_options(pfemSimulation PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(pfemSimulation PRIVATE /W4 /WX /wd4251)
endif()
