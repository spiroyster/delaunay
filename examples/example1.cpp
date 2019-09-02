#include "..\include\delaunay.hpp"

#include <sstream>
#include <fstream>
#include <random>

void writeOBJ(const std::string& filename, const std::list<Delaunay::Triangle>& triangles);

int main(int argc, char** argv)
{
	// Triangle...
	{
		Delaunay::Tessellator tessellator;

		std::vector<Delaunay::Vector2> vertices({
			Delaunay::Vector2(0, 0),
			Delaunay::Vector2(1.0, 0),
			Delaunay::Vector2(0.5, 1.0)
			});

		writeOBJ("triangle.obj", tessellator.triangulate(vertices));
	}
	

	// Sqaure...
	{
		Delaunay::Tessellator tessellator;

		std::vector<Delaunay::Vector2> vertices({
			Delaunay::Vector2(0, 0),
			Delaunay::Vector2(1.0, 0),
			Delaunay::Vector2(1.0, 1.0),
			Delaunay::Vector2(0, 1.0)
			});

		writeOBJ("square.obj", tessellator.triangulate(vertices));
	}

	// Ellipse...
	{
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
		
		writeOBJ("ellipse.obj", tessellator.triangulate(vertices));
	}

	// Constrained...
	{
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

		tessellator.triangulate(vertices);

		tessellator.addConstraint(std::pair<unsigned int, unsigned int>(0, segments / 2));

		writeOBJ("contrained.obj", tessellator.getTriangles());
	}

	// Constrained example 2...
	{
		Delaunay::Tessellator tessellator;

		std::vector<Delaunay::Vector2> vertices;
		std::vector<unsigned int> indexes;

		unsigned int n = 100;

		std::default_random_engine rng(n);
		std::uniform_real_distribution<double> rng_dist(0, 255);

		for (unsigned int p = 0; p < n; ++p)
		{
			vertices.push_back(Delaunay::Vector2(rng_dist(rng), rng_dist(rng)));
			indexes.push_back(p);
		}

		tessellator.triangulate(vertices);
		
		tessellator.addConstraint(std::pair<unsigned int, unsigned int>(0, n-1));
		
		

		writeOBJ("constrained2.obj", tessellator.getTriangles());
	}

	// Random points... 1000
	{
		Delaunay::Tessellator tessellator;

		std::vector<Delaunay::Vector2> vertices;

		unsigned int n = 1000;

		std::default_random_engine rng(n);
		std::uniform_real_distribution<double> rng_dist(0, 255);

		for (unsigned int p = 0; p < n; ++p)
			vertices.push_back(Delaunay::Vector2(rng_dist(rng), rng_dist(rng)));
		
		writeOBJ("random.obj", tessellator.triangulate(vertices));
	}


	// Append example...

	return 0;
}


void writeOBJ(const std::string& filename, const std::list<Delaunay::Triangle>& triangles)
{
	std::ostringstream oss;

	oss << "# delaunay output. https://github.com/spiroyster/delaunay\n";

	oss << "\n# vertices...\n";

	std::for_each(triangles.begin(), triangles.end(),
		[&oss](const Delaunay::Triangle& v)
	{
		oss << "v " << v.a->x << " " << v.a->y << " " << 0 << '\n';
		oss << "v " << v.b->x << " " << v.b->y << " " << 0 << '\n';
		oss << "v " << v.c->x << " " << v.c->y << " " << 0 << '\n';
	});

	oss << "\n# triangles...\n";

	unsigned int triangleIndex = 0;
	std::for_each(triangles.begin(), triangles.end(),
		[&oss, &triangleIndex](const Delaunay::Triangle&)
	{
		oss << "f " << ++triangleIndex << " " << ++triangleIndex << " " << ++triangleIndex << '\n';
	});

	std::ofstream file(filename);
	if (file)
		file << oss.str();
	else
		std::exception("Unable to save obj file.");
}