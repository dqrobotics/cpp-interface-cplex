CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(dqrobotics-interface-cplex)
set (CMAKE_CXX_STANDARD 11)

include_directory(include/dqrobotics/solvers)

################################################################
# INSTALL HEADERS
################################################################

INSTALL(FILES
    include/dqrobotics/solvers/DQ_CPLEXSolver.h
    DESTINATION "include/dqrobotics/solvers")
