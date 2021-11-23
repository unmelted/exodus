## 2021.09.16 
## created by Kelly (ekma@4dreplay.com)
##

cmake_minimum_required(VERSION 3.7)
#set(CMAKE_C_COMPILER "clang")
#set(CMAKE_CXX_COMPILER "/usr/bin/clang")

project(Exodus)
set (PLATFROM MAC) # available : WIN, ARM, LINUX
set(CMAKE_CXX_STANDARD 14)
find_package( OpenCV REQUIRED )

set(CMAKE_CXX_FLAGS "-Wno-invalid-source-encoding -Wno-deprecated-register -Wno-writable-strings -Wno-parentheses -O2 -g")
#set(LIB_INCLUDE_DIR "/usr/local/include")
add_definitions(-D_DEBUG)
add_definitions(-D_IMGDEBUG)

if(APPLE)
    add_definitions(-D_MAC_)
else ()
    add_definitions(-D_WIN_)
endif()

set(PRJ_PATH        ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB PRJ_SRC      ${PRJ_PATH}/src/*.cpp
                            ${PRJ_PATH}/src/common/*.cpp
                            ${PRJ_PATH}/exodus.cpp)

list(REMOVE_ITEM PRJ_SRC ${PRJ_PATH}/src/ESMConnector.cpp)

file(GLOB PRJ_INC      ${PRJ_PATH}/src/*.hpp
                            ${PRJ_PATH}/src/common/*.hpp
                            ${PRJ_PATH}/include/*.hpp
                            )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${PRJ_INC} ${OpenCV_INCLUDE_DIRS})
message(${OpenCV_INCLUDE_DIRS})
add_library(exodus SHARED ${PRJ_SRC})
target_link_libraries(exodus ${OpenCV_LIBS})

#add_executable(runtest ${PRJ_SRC})
#target_link_libraries(runtest)
