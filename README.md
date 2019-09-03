# About

Header only constrained delaunay tessellator, C++17 & STL only. 

This is originally a fork of https://github.com/Bl4ckb0ne/delaunay-triangulation (Bowyer-Watson) with a few bug fixes and added support for constrained triangulations. (so Kudos to Bl4ckb0ne for original implementation!).

# Usage

Simply include delaunay.hpp and you're off!

    #include "delaunay.hpp"

The Delaunay namespace provides a Tessellator class which you can use to create the tessellation of triangles given the input vertices. 

    Delaunay::Tessellator tessellator;

    std::vector<Delaunay::Vector2> vertices({
			Delaunay::Vector2(0, 0),
			Delaunay::Vector2(1.0, 0),
			Delaunay::Vector2(0.5, 1.0)
		});

    std::list<Delaunay::Triangle> result = tessellator.triangulate(vertices);

No stiener points are inserted so the number and ordering of vertices are preserved. The result is a list of triangles with pointers to the vertices contained in the tessellator class.

![Alt text](triangle.PNG?raw=true "Triangle")

Vertices must be unique and not duplicated otherwise this will cause duplicate triangles to be calculated. A function is provided for convinience for this.

    ValidateVerticesForTessellation(vertices);

however since this function sorts these vertices before uniquing them, the ordering may not be preservered for the original vertices. Do not use if this is important.

See examples/examples1.cpp for some more examples.

## Constrained tessellation

Constraints can be used to force edges in the resultant triangulation. First tessellate the vertices.

e.g This snippet tessellated vertices definiing a ellipse which is 2 units wide, 1 unit high, with 16 segments.

    Delaunay::Tessellator tessellator;

    unsigned int segments = 16;
    double width = 2.0;
    double height = 1.0;

    std::vector<Delaunay::Vector2> vertices;
    
    for (unsigned int d = 0; d < segments; ++d)
    {
        double angle = d * (2.0 * Delaunay::pi) / segments;
        vertices.push_back(Delaunay::Vector2(width * cos(angle), height * sin(angle)));
    }

    std::list<Delaunay::Triangle> result = tessellator.triangulate(vertices);

![Alt text](ellipse.PNG?raw=true "Ellipse")

Constraines are added as vertex index pairs. Constraints are added individually, a single constraint at a time. If the edge defined by the two vertex indexes is already present in the tessellation, it is ignored.

    tessellator.addConstraint(std::pair<unsigned int, unsigned int>(0, segments / 2));

Then call getTriangles() to get the resultant triangles from the constrained tessellation.

![Alt text](constrained.PNG?raw=true "Constrained")

When adding multiple constraints, edges defining the constraints should ideally not intersect each other, however since they are added individually, one constraining edge will overwrite the other.

e.g Adding another constraint to the above tessellation which crosses the original constraint will result in the second constraint given priority, effectively loosing the previous one.

    tessellator.addConstraint(std::pair<unsigned int, unsigned int>(segments / 4, 3 * (segments / 4)));

![Alt text](constrained2.PNG?raw=true "Constrained2")

See examples1.cpp for the above code in all its glory.

# Bugs

Please submit bug data (if applicable) as .csv. A CSV output routiune is provided for convinience.

    Delaunay::csv(...);

This does not create a csv file, only the syntax for a csv file, so caller is reposonsible for writing a file.

# Future

* Potentially remove std::optional to allow compatiblilty with earlier C++ standards. std::optional is a C++17 addition.
* Add tests
* Doxygenate

![Alt text](footer.PNG?raw=true "Footer")