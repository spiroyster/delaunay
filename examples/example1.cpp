#include "..\include\delaunay.hpp"

#include <sstream>
#include <fstream>
#include <random>

void writeOBJ(const std::string& filename, const delaunay::tessellator& t);

int main(int argc, char** argv)
{
	// Triangle...
	{
		std::vector<delaunay::vector2> vertices({
			delaunay::vector2(0, 0),
			delaunay::vector2(1.0, 0),
			delaunay::vector2(0.5, 1.0)
			});

		delaunay::tessellator t(vertices);
		writeOBJ("triangle.obj", t);
	}

	// Sqaure...
	{
		std::vector<delaunay::vector2> vertices({
			delaunay::vector2(0, 0),
			delaunay::vector2(1.0, 0),
			delaunay::vector2(1.0, 1.0),
			delaunay::vector2(0, 1.0)
			});

		delaunay::tessellator t(vertices);
		writeOBJ("square.obj", t);
	}


	// Ellipse...
	{
		unsigned int segments = 16;
		double width = 2.0;
		double height = 1.0;

		std::vector<delaunay::vector2> vertices;
		
		for (unsigned int d = 0; d < segments; ++d)
		{
			double angle = d * (2.0 * delaunay::pi) / segments;
			vertices.push_back(delaunay::vector2(width * cos(angle), height * sin(angle)));
		}
		
		delaunay::tessellator t(vertices);
		writeOBJ("ellipse.obj", t);
	}

	// Constrained...
	{
		unsigned int segments = 16;
		double width = 2.0 /2.0;
		double height = 1.0 /2.0;
		double twoPI = 2.0 * delaunay::pi;

		std::vector<delaunay::vector2> vertices;

		for (unsigned int d = 0; d < segments; ++d)
		{
			double angle = d * twoPI / segments;
			vertices.push_back(delaunay::vector2(width * cos(angle), height * sin(angle)));
		}

		delaunay::tessellator t(vertices);

		t.addConstraint(delaunay::edge{ 0, segments / 2 });

		writeOBJ("constrained.obj", t);

		t.addConstraint(delaunay::edge{ segments / 4, 3 * (segments / 4) });

		writeOBJ("constrained2.obj", t);
	}

	// Constrained with steiner...
	{
		unsigned int segments = 16;
		double width = 2.0 / 2.0;
		double height = 1.0 / 2.0;
		double twoPI = 2.0 * delaunay::pi;

		std::vector<delaunay::vector2> vertices;

		for (unsigned int d = 0; d < segments; ++d)
		{
			double angle = d * twoPI / segments;
			vertices.push_back(delaunay::vector2(width * cos(angle), height * sin(angle)));
		}

		// insert a steiner point in the centre of the ellipse...
		vertices.push_back(delaunay::vector2(0, 0));
		unsigned int steinerIndex = vertices.size() - 1;

		delaunay::tessellator t(vertices);

		writeOBJ("constrained3.obj", t);

		t.addConstraint(delaunay::edge{ 0, steinerIndex });
		t.addConstraint(delaunay::edge{ steinerIndex, segments / 2 });
		
		// N.B Due to the default tessellation done, these constraints should already be present as edges, stated here for completness...
		t.addConstraint(delaunay::edge{ segments / 4, steinerIndex });
		t.addConstraint(delaunay::edge{ steinerIndex, 3 * (segments / 4) });

		writeOBJ("constrained4.obj", t);
	}

	// Random points... 1000
	{
		std::vector<delaunay::vector2> vertices;

		unsigned int n = 1000;

		std::default_random_engine rng(n);
		std::uniform_real_distribution<double> rng_dist(0, 255);

		for (unsigned int p = 0; p < n; ++p)
			vertices.push_back(delaunay::vector2(rng_dist(rng), rng_dist(rng)));
		
		delaunay::tessellator t(vertices);
		writeOBJ("random.obj", t);
	}

	return 0;
}


void writeOBJ(const std::string& filename, const delaunay::tessellator& tessellator)
{
	std::ostringstream oss;

	oss << "# delaunay output. https://github.com/spiroyster/delaunay\n";

	oss << "\n# vertices...\n";

	std::for_each(tessellator.getVertices().begin(), tessellator.getVertices().end(),
		[&oss, &tessellator](const delaunay::vector2& v)
	{
		oss << "v " << v.x_ << " " << v.y_ << " " << 0 << '\n';
	});

	oss << "\n# triangles...\n";

	unsigned int triangleIndex = 0;
	std::for_each(tessellator.getTriangles().begin(), tessellator.getTriangles().end(),
		[&oss](const delaunay::triangle& t)
	{
		oss << "f " << t.a_+1 << " " << t.b_+1 << " " << t.c_+1 << '\n';
	});

	std::ofstream file(filename);
	if (file)
		file << oss.str();
	else
		std::exception("Unable to save obj file.");
}