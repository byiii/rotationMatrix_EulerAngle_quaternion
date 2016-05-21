cmake_minimum_required(VERSION 2.8)
project(rotationMatrix)

set(Eigen_INCLUDE_DIR D:/Development/Eigen/include/eigen3)
INCLUDE_DIRECTORIES(${Eigen_INCLUDE_DIR})

aux_source_directory(. SRC_LIST)

add_executable(${PROJECT_NAME} ${SRC_LIST})
