find_path(EIGEN_INCLUDE_DIRS NAMES "Eigen")
if(NOT EIGEN_INCLUDE_DIRS)
    message(FATAL_ERROR "Eigen not found!")
else()
	message(STATUS "Found Eigen: " ${EIGEN_INCLUDE_DIRS})
endif()

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "Found CGAL: " ${CGAL_DIR})
endif()

if(USE_TBB_CGAL)
    find_package(TBB)
    if(TBB_FOUND)
        message(STATUS "Found TBB: " ${TBB_LIBRARIES})
        include(CGAL_TBB_support)
    else()
        message(FATAL_ERROR "TBB not found!")
    endif()
endif()

file(GLOB_RECURSE LIBMESH_SRCS
     "${PROJECT_SOURCE_DIR}/srcs/mesh/*.h"
     "${PROJECT_SOURCE_DIR}/srcs/mesh/*.inl"
     "${PROJECT_SOURCE_DIR}/srcs/mesh/*.hpp"
     "${PROJECT_SOURCE_DIR}/srcs/mesh/*.cpp")

add_library(pfemMesh SHARED ${LIBMESH_SRCS})
target_include_directories(pfemMesh SYSTEM
                           PRIVATE ${EIGEN_INCLUDE_DIRS}
                           PRIVATE ${GMSH_INCLUDE_DIRS})
target_link_libraries(pfemMesh PRIVATE OpenMP::OpenMP_CXX  CGAL::CGAL ${GMSH_LIBRARIES})
if(USE_TBB_CGAL)
    target_link_libraries(pfemMesh PRIVATE CGAL::TBB_support)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    target_compile_options(pfemMesh PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    target_compile_options(pfemMesh PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(pfemMesh PRIVATE /W4 /WX /wd4251)
endif()
