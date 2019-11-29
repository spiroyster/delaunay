#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "objio.hpp"

#include "..\include\delaunay.hpp"

namespace
{
	std::vector<delaunay::vector2> readCSV(const std::string& filename)
	{
		std::vector<delaunay::vector2> result;

		std::ifstream syntax(filename);
		if (!syntax)
			throw std::exception(std::string("Unable to read file " + filename).c_str());

		std::string line;
		unsigned int lineNum = 1;

		while (std::getline(syntax, line))
		{
			std::size_t comma = line.find(", ");

			if (comma == std::string::npos)
				throw std::exception(std::string(filename + "(" + std::to_string(lineNum) + "):" + line).c_str());

			result.push_back(delaunay::vector2(std::stod(line.substr(0, comma)), std::stod(line.substr(comma + 1))));
			++lineNum;
		}

		return result;
	}

	std::vector<delaunay::triangle> readOBJ(const std::string& filename)
	{
		std::shared_ptr<objio::mesh> obj = objio::readFile(filename);
		std::vector<delaunay::triangle> result;
		result.reserve(obj->faces_.size());

		std::for_each(obj->faces_.begin(), obj->faces_.end(), 
			[&result](const objio::face& f) 
			{
				result.push_back(delaunay::triangle{ static_cast<unsigned int>(f.vertices_[0].p_), static_cast<unsigned int>(f.vertices_[1].p_), static_cast<unsigned int>(f.vertices_[2].p_) });
			});

		return result;
	}

	bool approxEqual(const double& x, const double& y, int ulp = 2)
	{
		// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
		return std::abs(x - y) <= std::numeric_limits<double>::epsilon() * std::abs(x + y) * static_cast<double>(ulp) || std::abs(x - y) < (std::numeric_limits<double>::min)();
	}	
	bool approxEqual(const delaunay::vector2& a, const delaunay::vector2& b)
	{
		return approxEqual(a.x_, b.x_) && approxEqual(a.y_, b.y_);
	}
	bool equal(const delaunay::triangle& t, const delaunay::triangle& tt)
	{
		return (t.a_ == tt.a_ || t.a_ == tt.b_ || t.a_ == tt.c_) && (t.b_ == tt.a_ || t.b_ == tt.b_ || t.b_ == tt.c_) && (t.c_ == tt.a_ || t.c_ == tt.b_ || t.c_ == tt.c_);
	}
		

	void check(const std::string& csvFilename, const std::string& objFilename)
	{
		// Load in the csv...
		std::vector<delaunay::vector2> vertices = readCSV(csvFilename);

		// Load in the obj...
		std::vector<delaunay::triangle> triangles = readOBJ(objFilename);

		// perform the tessellation in accordance with the csv vertices...
		delaunay::tessellator tessellator(vertices);

		// compare the vertices...
		REQUIRE(tessellator.getVertices().size() == vertices.size());
		for (unsigned int i = 0; i < vertices.size(); ++i)
		{
			const delaunay::vector2& v = vertices[i];
			REQUIRE(std::find_if(tessellator.getVertices().begin(), tessellator.getVertices().end(), [&v](const delaunay::vector2& vv) { return approxEqual(v, vv); }) != tessellator.getVertices().end());
		}

		// compare to the obj triangles
		REQUIRE(tessellator.getTriangles().size() == triangles.size());
		for (unsigned int i = 0; i < triangles.size(); ++i)
		{
			const delaunay::triangle& t = triangles[i];
			REQUIRE(std::find_if(tessellator.getTriangles().begin(), tessellator.getTriangles().end(), [&t](const delaunay::triangle& tt) { return equal(t, tt); }) != tessellator.getTriangles().end());
		}

	}

	void check(const std::string& csvFilename, const std::string& edgesFilename, const std::string& objFilename)
	{
		// Load in the csv...
		std::vector<delaunay::vector2> vertices = readCSV(csvFilename);

		// Load in the edges csv...
		std::vector<delaunay::vector2> edges = readCSV(edgesFilename);

		// Load in the obj...
		std::vector<delaunay::triangle> triangles = readOBJ(objFilename);

		// perform the tessellation in accordance with the csv vertices...
		delaunay::tessellator tessellator(vertices);

		std::for_each(edges.begin(), edges.end(), 
			[&tessellator](const delaunay::vector2& edge) 
			{
				tessellator.addConstraint(delaunay::edge{ static_cast<unsigned int>(edge.x_), static_cast<unsigned int>(edge.y_) });
			});

		// compare the vertices...
		REQUIRE(tessellator.getVertices().size() == vertices.size());
		for (unsigned int i = 0; i < vertices.size(); ++i)
		{
			const delaunay::vector2& v = vertices[i];
			REQUIRE(std::find_if(tessellator.getVertices().begin(), tessellator.getVertices().end(), [&v](const delaunay::vector2& vv) { return approxEqual(v, vv); }) != tessellator.getVertices().end());
		}

		// compare to the obj triangles
		REQUIRE(tessellator.getTriangles().size() == triangles.size());
		for (unsigned int i = 0; i < triangles.size(); ++i)
		{
			const delaunay::triangle& t = triangles[i];
			REQUIRE(std::find_if(tessellator.getTriangles().begin(), tessellator.getTriangles().end(), [&t](const delaunay::triangle& tt) { return equal(t, tt); }) != tessellator.getTriangles().end());
		}

	}
}

TEST_CASE("Tessellate triangle", "[delaunay_tessellate_triangle]") { check("data/triangle.csv", "data/triangle.obj"); }
TEST_CASE("Tessellate square", "[delaunay_tessellate_sqaure]") { check("data/square.csv", "data/square.obj"); }
TEST_CASE("Tessellate ellipse", "[delaunay_tessellate_ellipse]") { check("data/ellipse.csv", "data/ellipse.obj"); }
TEST_CASE("Tessellate ellipse with constraint", "[delaunay_tessellate_ellipse_constaint]") { check("data/ellipse.csv", "data/ellipse_constraint.csv", "data/ellipse_constraint.obj"); }
TEST_CASE("Tessellate ellipse with constraint overwrite", "[delaunay_tessellate_ellipse_constaint_overwrite]") { check("data/ellipse.csv", "data/ellipse_constraint_overwrite.csv", "data/ellipse_constraint_overwrite.obj"); }

