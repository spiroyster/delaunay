#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "..\include\delaunay.hpp"

bool approxEqual(const double& x, const double& y, int ulp = 2)
{
	// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
	return std::abs(x - y) <= std::numeric_limits<double>::epsilon() * std::abs(x + y) * static_cast<double>(ulp) || std::abs(x - y) < (std::numeric_limits<double>::min)();
}

bool approxEqual(const delaunay::vector2& a, const delaunay::vector2& b)
{
	return approxEqual(a.x_, b.x_) && approxEqual(a.y_, b.y_);
}

bool equal(const delaunay::edge& e, const delaunay::edge& ee)
{
	return (e.i_ == ee.i_ && e.j_ == ee.j_) || (e.i_ == ee.j_ && e.j_ == ee.i_);
}

bool equal(const delaunay::triangle& t, const delaunay::triangle& tt)
{
	return (t.a_ == tt.a_ || t.a_ == tt.b_ || t.a_ == tt.c_) && (t.b_ == tt.a_ || t.b_ == tt.b_ || t.b_ == tt.c_) && (t.c_ == tt.a_ || t.c_ == tt.b_ || t.c_ == tt.c_);
}

void check(const delaunay::tessellator& tessellator, const std::vector<delaunay::vector2>& vertices, const std::vector<delaunay::edge>& edges, const std::vector<delaunay::triangle>& triangles)
{
	// Check the vertices...
	REQUIRE(tessellator.getVertices().size() == vertices.size());
	for (unsigned int i = 0; i < vertices.size(); ++i)
	{
		const delaunay::vector2& v = vertices[i];
		REQUIRE(std::find_if(tessellator.getVertices().begin(), tessellator.getVertices().end(), [&v](const delaunay::vector2& vv) { return approxEqual(v, vv); }) != tessellator.getVertices().end());
	}
		
	REQUIRE(tessellator.getEdges().size() == edges.size());
	for (unsigned int i = 0; i < edges.size(); ++i)
	{
		const delaunay::edge& e = edges[i];
		REQUIRE(std::find_if(tessellator.getEdges().begin(), tessellator.getEdges().end(), [&e](const delaunay::edge& ee) { return equal(e, ee); }) != tessellator.getEdges().end());
	}

	REQUIRE(tessellator.getTriangles().size() == triangles.size());
	for (unsigned int i = 0; i < triangles.size(); ++i)
	{
		const delaunay::triangle& t = triangles[i];
		REQUIRE(std::find_if(tessellator.getTriangles().begin(), tessellator.getTriangles().end(), [&t](const delaunay::triangle& tt) { return equal(t, tt); }) != tessellator.getTriangles().end());
	}

}

TEST_CASE("Tessellate triangle", "[delaunay_tessellate_triangle]")
{
	//Delaunay::Tessellator tessellator;

	// Vertices...
	std::vector<delaunay::vector2> vertices({
		delaunay::vector2(0, 0),
		delaunay::vector2(1.0, 0),
		delaunay::vector2(0.5, 1.0)
		});

	// Edges...
	std::vector<delaunay::edge> edges({
		delaunay::edge{0, 1},
		delaunay::edge{1, 2},
		delaunay::edge{2, 0}
		});
	
	// Triangles...
	std::vector<delaunay::triangle> triangles({
		delaunay::triangle{0, 1, 2}
		});

	delaunay::tessellator t(vertices);

	check(t, vertices, edges, triangles);
}

TEST_CASE("Tessellate square", "[delaunay_tessellate_square]")
{
	// Vertices...
	std::vector<delaunay::vector2> vertices({
			delaunay::vector2(0, 0),
			delaunay::vector2(1.0, 0),
			delaunay::vector2(1.0, 1.0),
			delaunay::vector2(0, 1.0)
		});

	// Edges...
	std::vector<delaunay::edge> edges({
		delaunay::edge{0, 1},
		delaunay::edge{1, 2},
		delaunay::edge{2, 3},
		delaunay::edge{3, 0},
		delaunay::edge{1, 3}
		});

	// Triangles...
	std::vector<delaunay::triangle> triangles({
		delaunay::triangle{0, 1, 3},
		delaunay::triangle{1, 2, 3}
		});

	delaunay::tessellator t(vertices);

	check(t, vertices, edges, triangles);
}