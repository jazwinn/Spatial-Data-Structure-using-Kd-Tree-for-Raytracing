# k-d Tree (KDTree) Demo

Demonstration of a k-d Tree (KDTree) for spatial partitioning to accelerate ray tracing operations, including fast intersection tests and visibility determination.

## Demo

![1758038139840](https://github.com/user-attachments/assets/0237780a-c774-4383-a3e0-8df24c10fd5e)

## Overview

A **k-d Tree** is a binary tree structure used for organizing points or objects in a k-dimensional space. Each node represents a **split** along one dimension and partitions the space into two subspaces. KDTrees are widely used in ray tracing for:

* **Fast ray-object intersection tests**.
* **Efficient scene traversal**.
* **Visibility determination and shading acceleration**.

The main idea is to **reduce the number of intersection tests** by recursively partitioning the space and skipping large subspaces that a ray does not intersect.

### Prerequisites
Before you begin, make sure you have the following installed:
- **Cmake**
- **Cmake with visual studio**
- **vcpkg**

## Install Packages via vcpkg

**Windows (x64):**

```powershell
.\vcpkg\vcpkg install glm:x64-windows
.\vcpkg\vcpkg install glad:x64-windows
.\vcpkg\vcpkg install glfw3:x64-windows
.\vcpkg\vcpkg install lodepng:x64-windows
.\vcpkg\vcpkg install gtest:x64-windows
.\vcpkg\vcpkg install fmt:x64-windows
.\vcpkg\vcpkg install imgui[glfw-binding,opengl3-binding]:x64-windows
.\vcpkg\vcpkg install imguizmo:x64-windows
.\vcpkg\vcpkg install stb:x64-windows
```
  
**Linux (x64):**
```bash
   ./vcpkg/vcpkg install glm
   ./vcpkg/vcpkg install glad
   ./vcpkg/vcpkg install glfw3
   ./vcpkg/vcpkg install lodepng
   ./vcpkg/vcpkg install gtest
   ./vcpkg/vcpkg install fmt
   ./vcpkg/vcpkg install imgui[glfw-binding,opengl3-binding]
   ./vcpkg/vcpkg install imguizmo
   ./vcpkg/vcpkg install stb
```

## Steps

1. **Clone the Repository**

```bash
git clone https://github.com/<your-username>/Spatial-Data-Structure-using-Bounding-Volume-Hierarchy-BVH-.git
cd Spatial-Data-Structure-using-Bounding-Volume-Hierarchy-BVH-
```

2. **Create a Build Directory**

```bash
mkdir build
cd build
```

3. **Configure the Project with CMake**

```bash
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake
```

4. **Build the Project**

* **Windows (Visual Studio)**

```powershell
cmake --build . --config Release
```

* **Linux**

```bash
cmake --build .
```

5. **Run the Executable**

```bash
./cs350-demo.exe
```


## Side Note
MODEL ASSETS REMOVED FOR COPYRIGHT REASON, PM FOR MORE DETAIL
You can see visual representation of a bounding volume hierachy under "BVH project\Visual Graph"


Graphs .dot or .txt : "\cs350su25_jazwinn.ng_a3\Submission Resource\dot&txt"
Graphs .png : "\cs350su25_jazwinn.ng_a3\Submission Resource\png"
Graphs .svg : "\cs350su25_jazwinn.ng_a3\Submission Resource\svg"

Switch Camera/View Frustrum Camera: 
1. Imgui Panel "BVH"
2. Header "Camera"
3. Select "Switch to Frustrum" to true to control the frustum, 
   switch back (false) to gain control of the other camera and see the frustum.
   Select "Draw Frustrum" to turn drawing frustrum on and off.

   
