cmake_minimum_required(VERSION 3.8)
project(cs350-demo)

# Demo
add_executable(${PROJECT_NAME} main.cpp ./ImGui.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC cs350-lib)

# ImGui
find_package(imgui CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC imgui::imgui)

# ImGuizmo
find_package(imguizmo CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE imguizmo::imguizmo)
