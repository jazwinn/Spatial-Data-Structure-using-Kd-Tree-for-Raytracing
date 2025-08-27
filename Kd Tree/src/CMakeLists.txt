cmake_minimum_required(VERSION 3.11)
project(cs350-lib)

# Engine library
add_library(
        ${PROJECT_NAME}
        Camera.cpp
        Camera.hpp
        Window.hpp
        Window.cpp
        OpenGl.hpp
        OpenGl.cpp
        Utils.hpp
        Utils.cpp
        Shader.hpp
        Shader.cpp
        Primitive.hpp
        Primitive.cpp
        CS350Loader.hpp
        CS350Loader.cpp
        Logging.hpp
        Logging.cpp
        PRNG.cpp
        PRNG.h
        Math.hpp
        DebugRenderer.hpp
        DebugRenderer.cpp
        KdTree.hpp
        Shapes.hpp
        Shapes.cpp
        Stats.hpp )
target_include_directories(${PROJECT_NAME} PUBLIC .)

# GLM
find_package(lodepng CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC lodepng)

# GLM
find_package(glm CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC glm::glm)

# GLAD
find_package(glad CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC glad::glad)

# GLFW3
find_package(glfw3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC glfw)

# fmt
find_package(fmt CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC fmt::fmt)
