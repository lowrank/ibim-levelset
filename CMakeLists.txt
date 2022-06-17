cmake_minimum_required(VERSION 2.8)
project(levelset)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -g -Ofast -DRUN_OMP -march=native -std=c++11 -fopenmp ")

find_package(GLUT)

if (GLUT_FOUND)
    set(SOURCE_FILES main.cpp levelset.cpp levelset.h utils.h molecule.h ls_point.h view.h view.cpp linalg.h blas_wrapper.h lapacke_wrapper.h gmres.h bbfmm.h bbfmm.cpp linalg.cpp blas_wrapper.cpp lapacke_wrapper.cpp electric.cpp electric_correction.cpp gmres.cpp Config.cpp Config.h)
    add_executable(levelset ${SOURCE_FILES})
    target_link_libraries(levelset m GL GLU glut gfortran ${CMAKE_SOURCE_DIR}/contrib/lib/libopenblas.a)
else ()
    set(SOURCE_FILES main.cpp levelset.cpp levelset.h utils.h molecule.h ls_point.h linalg.h blas_wrapper.h lapacke_wrapper.h gmres.h bbfmm.h bbfmm.cpp blas_wrapper.cpp lapacke_wrapper.cpp linalg.cpp electric.cpp electric_correction.cpp gmres.cpp Config.cpp Config.h)
    add_executable(levelset ${SOURCE_FILES})
    target_link_libraries(levelset m gfortran ${CMAKE_SOURCE_DIR}/contrib/lib/libopenblas.a)
endif ()

