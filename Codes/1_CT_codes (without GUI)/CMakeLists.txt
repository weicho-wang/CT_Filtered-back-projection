cmake_minimum_required(VERSION 3.16)
project(CT_Filtered_Back_Projection LANGUAGES CXX)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(OpenCV_DIR "C:\\opencv\\opencv\\build\\x64\\vc15\\lib")

# 查找 OpenCV
find_package(OpenCV REQUIRED)

# 收集所有源文件
set(SOURCES
    CT_main.cpp
    basic_display_functions.cpp
    basic_functions.cpp
    basic_radon_functions.cpp
    filter_parallel_line.cpp
    radon_fan_beam.cpp
    radon_parallel_line.cpp
    shepplogan.cpp
)

# 创建可执行文件
add_executable(CT_Reconstruction ${SOURCES})

# 链接 OpenCV
target_link_libraries(CT_Reconstruction ${OpenCV_LIBS})

# 包含头文件目录（当前目录）
target_include_directories(CT_Reconstruction PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

# 如果是 MSVC，设置一些选项
if(MSVC)
    target_compile_options(CT_Reconstruction PRIVATE /W3 /sdl)
    set_property(TARGET CT_Reconstruction PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

set_property(TARGET CT_Reconstruction PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})