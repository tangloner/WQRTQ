cmake_minimum_required(VERSION 3.12)

project(WQRTQ)

set(CMAKE_CXX_STANDARD 14)

include_directories(QuadProg/include)
include_directories(RTree/include)


file(GLOB SOURCES "*.h")
file(GLOB LIBRTREE "RTree/*.cpp")
file(GLOB LIBQUADPROG "QuadProg/*.cc")

add_executable(WQRTQ main.cpp ${SOURCES} ${LIBQUADPROG} ${LIBRTREE} minheap.h)
