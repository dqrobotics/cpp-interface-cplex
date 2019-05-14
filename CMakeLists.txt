CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(dqrobotics-interface-cplex)
set (CMAKE_CXX_STANDARD 11)

#ADD_DEFINITIONS(-g -O2 -Wall)
FIND_PACKAGE(Eigen3 REQUIRED)

add_definitions(-DIL_STD)

# CPLEX
SET(CPLEX_DIR "/opt/ibm/ILOG/CPLEX_Studio1271")
INCLUDE_DIRECTORIES(${CPLEX_DIR}/concert/include)
INCLUDE_DIRECTORIES(${CPLEX_DIR}/cplex/include)

# Eigen3
INCLUDE_DIRECTORIES(EIGEN3_INCLUDE_DIR)
# This project
INCLUDE_DIRECTORIES(include)

################################################################
# DEFINE AND INSTALL LIBRARY AND INCLUDE FOLDER
################################################################

ADD_LIBRARY(${PROJECT_NAME} SHARED 
    src/cplex_interface.cpp
    )

SET_TARGET_PROPERTIES(${PROJECT_NAME} 
    PROPERTIES PUBLIC_HEADER
    include/dqrobotics/interfaces/cplex_interface.h
    )

INSTALL(TARGETS ${PROJECT_NAME} 
    LIBRARY DESTINATION "lib"
    PUBLIC_HEADER DESTINATION "include/dqrobotics/interfaces"
    PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

################################################################
# INSTALL SOURCE FILES (So that the debugger can find them)
################################################################

# base folder
INSTALL(FILES 
    src/cplex_interface.cpp
    DESTINATION "src/dqrobotics/interfaces")
