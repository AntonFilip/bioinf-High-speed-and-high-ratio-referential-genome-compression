cmake_minimum_required(VERSION 3.15)
project(bioinf_High_speed_and_high_ratio_referential_genome_compression)

set(CMAKE_CXX_STANDARD 14)

add_executable(hirgc_comp hirgc_comp.cpp)
add_executable(hirgc_decomp hirgc_decomp.cpp)