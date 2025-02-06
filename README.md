# graphics-engine
A 2D graphics engine featuring clipping, polygon and bezier curve support, circle rendering, path manipulation, eight blend modes, and three custom shaders.

## Features
*   **Clipping:** Supports clipping of drawn elements to specified rectangular bounds.
*   **Polygons:** Renders filled convex polygons with specified colors and blend modes.  Handles arbitrary numbers of vertices.
*   **Bezier Curves:** Draws quadratic and cubic Bezier curves with customizable precision.
*   **Circles:** Efficiently renders circles using quadratic curves.
*   **Paths:** Uses a flexible path API (`GPath`) allowing for complex shapes constructed from lines, curves, and rectangles.
*   **Blend Modes:** Implements eight different blend modes (`GBlendMode`) for compositing layers.
*   **Shaders:** Includes three custom shaders:
    *   `MyShader`:  A bitmap shader supporting clamping, repeating, and mirroring tile modes.
    *   `MyVoronoiShader`: Creates a Voronoi diagram effect based on input points and colors.
    *   `MyLinearPosGradientShader`: Generates linear gradients with colors positioned along a line based on specified percentages.
*   **Mesh Rendering:**  Draws arbitrary triangle meshes with optional per-vertex colors and texture coordinates.

## Installation
1.  **Clone the repository:** `git clone https://github.com/aden-cle/graphics-engine.git`

## Usage
1. Edit the `GDrawSomething` function in `MyCanvas.cpp` by performing drawing operations (`drawPath`, `drawMesh`, etc.) on the provided canvas.

2. **Compile:** Use the provided Makefile: `make`.

3. **Run:** `./image -e expected -v` to generate your drawing.

**Example:**

```c++
// Make the background red
GColor red {1.0, 0.0, 0.0, 1.0};
canvas->clear(red);

// Create a circle with center at (100, 100) and radius 40
GPath path;
path = GPath();
path.addCircle({100, 100}, 40, GPath::kCCW_Direction);

// Make the circle green and draw it on the canvas
GColor green {0.0, 1.0, 0.0, 1.0};
GPaint green_paint(green);
canvas->drawPath(path, green_paint); 

return "My title here"
```

## Configuration
No specific configuration is required beyond compiling the code using the provided Makefile.

## API Documentation
Refer to the header files (`*.h`) for detailed function specifications and parameters.
