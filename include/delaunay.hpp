/*
	MIT License

	Copyright (c) 2019 Cordell Barron

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.

*/
#ifndef DELAUNAY_TRIANGULATION_HPP
#define DELAUNAY_TRIANGULATION_HPP

#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include <optional>
#include <sstream>
#include <assert.h>

namespace delaunay
{
	//const static double pi = 3.14159265358979323846;

	struct vector2 { vector2(double x, double y) : x_(x), y_(y) {} double x_, y_; };
	struct edge { unsigned int i_, j_; };
	struct triangle { unsigned int a_, b_, c_; };

	class tessellator
	{
		struct workingVertex
		{
			workingVertex(const vector2* v, unsigned int index)
				:	v_(v), index_(index)
			{
			}

			double dist(const vector2 &v) const
			{
				const double dx = v_->x_ - v.x_;
				const double dy = v_->y_ - v.y_;
				return dx * dx + dy * dy;
			}

			const vector2* v_;
			unsigned int index_;
		};

		struct workingTriangle
		{
			workingTriangle(const workingVertex* a, const workingVertex* b, const workingVertex* c)
				: a_(a), b_(b), c_(c), valid_(true) {}

			bool containsVertex(unsigned int index)
			{
				return a_->index_ == index || b_->index_ == index || c_->index_ == index;
			}

			const workingVertex* a_;
			const workingVertex* b_;
			const workingVertex* c_;
			bool valid_;

			bool isWithinCircumcircle(const workingVertex& v) const
			{
				const double ab = a_->v_->x_ * a_->v_->x_ + a_->v_->y_ * a_->v_->y_;
				const double cd = b_->v_->x_ * b_->v_->x_ + b_->v_->y_ * b_->v_->y_;
				const double ef = c_->v_->x_ * c_->v_->x_ + c_->v_->y_ * c_->v_->y_;

				const double ax = a_->v_->x_;
				const double ay = a_->v_->y_;
				const double bx = b_->v_->x_;
				const double by = b_->v_->y_;
				const double cx = c_->v_->x_;
				const double cy = c_->v_->y_;

				const double circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
				const double circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));

				vector2 circum(circum_x * 0.5, circum_y * 0.5);
				return  v.dist(circum) <= a_->dist(circum);
			}
		};

		struct workingEdge
		{
			workingEdge(const workingVertex* i, const workingVertex* j)
				: i_(i), j_(j), valid_(true)
			{
			}

			bool operator==(const workingEdge& rhs)
			{
				return (rhs.i_ == i_ && rhs.j_ == j_) || (rhs.i_ == j_ && rhs.j_ == i_);
			}

			const workingVertex* i_;
			const workingVertex* j_;
			bool valid_;
		};
		

	public:
		tessellator(const std::vector<vector2>& points)
			: vertices_(points)
		{
			std::vector<unsigned int> allVertexIndexes(points.size(), 0);
			for (unsigned int i = 0; i < points.size(); ++i)
				allVertexIndexes[i] = i;
			
			tessellate(allVertexIndexes);
		}

		void addConstraint(const edge& e)
		{
			// check if we have this constraint edge already...

			// find out which edges this constraint straddles...

			// retessellate both the right and left side
		}

		const std::list<triangle>& getTriangles() const { return triangles_; }
		const std::vector<vector2>& getVertices() const { return vertices_; }
		const std::list<edge>& getEdges() const { return edges_; }

	private:
		void tessellate(const std::vector<unsigned int>& vertexIndexes)
		{
			assert(!vertexIndexes.empty());
			assert(!vertices_.empty());

			// create our sublist that we tessellate to....
			std::vector<workingVertex> workingVertices(vertexIndexes.size() + 3, workingVertex(0, 0));
			for (unsigned int i = 0; i < vertexIndexes.size(); ++i)
				workingVertices[i] = workingVertex(&vertices_[vertexIndexes[i]], vertexIndexes[i]);

			// determine the super triangle...
			double minX = workingVertices[0].v_->x_;
			double minY = workingVertices[0].v_->y_;
			double maxX = minX;
			double maxY = minY;

			for (std::size_t i = 0; i < vertexIndexes.size(); ++i)
			{
				if (workingVertices[i].v_->x_ < minX) minX = workingVertices[i].v_->x_;
				if (workingVertices[i].v_->y_ < minY) minY = workingVertices[i].v_->y_;
				if (workingVertices[i].v_->x_ > maxX) maxX = workingVertices[i].v_->x_;
				if (workingVertices[i].v_->y_ > maxY) maxY = workingVertices[i].v_->y_;
			}

			double dx = maxX - minX;
			double dy = maxY - minY;
			double deltaMax = (std::max)(dx, dy);
			double midx = (minX + maxX) * 0.5;
			double midy = (minY + maxY) * 0.5;

			// Add our three super triangle vertices...
			unsigned int p1Index = vertexIndexes.size(), p2Index = vertexIndexes.size() + 1, p3Index = vertexIndexes.size() + 2;

			vector2 p1(midx - 20.0 * deltaMax, midy - deltaMax);
			vector2 p2(midx, midy + 20.0 * deltaMax);
			vector2 p3(midx + 20.0 * deltaMax, midy - deltaMax);

			workingVertices[vertexIndexes.size()] = workingVertex(&p1, p1Index);
			workingVertices[vertexIndexes.size()+1] = workingVertex(&p2, p2Index);
			workingVertices[vertexIndexes.size()+2] = workingVertex(&p3, p3Index);

			// Create a list of triangles, and add the supertriangle in it
			std::list<workingTriangle> workingTriangles;
			workingTriangles.push_back(workingTriangle(&workingVertices[p1Index], &workingVertices[p2Index], &workingVertices[p3Index]));

			std::for_each(workingVertices.begin(), workingVertices.end(),
				[&workingTriangles](workingVertex& v) 
				{
					std::vector<workingEdge> polygon;

					std::for_each(workingTriangles.begin(), workingTriangles.end(), [&v, &polygon](workingTriangle& t) 
						{
							if (t.isWithinCircumcircle(v))
							{
								t.valid_ = false;
								polygon.push_back(workingEdge(t.a_, t.b_));
								polygon.push_back(workingEdge(t.b_, t.c_));
								polygon.push_back(workingEdge(t.c_, t.a_));
							}
						});

					workingTriangles.erase(std::remove_if(workingTriangles.begin(), workingTriangles.end(), [](workingTriangle& t) { return !t.valid_; }), workingTriangles.end());

					for (auto e1 = polygon.begin(); e1 != polygon.end(); ++e1)
					{
						for (auto e2 = e1 + 1; e2 != polygon.end(); ++e2)
						{
							if (*e1 == *e2)
							{
								e1->valid_ = false;
								e2->valid_ = false;
							}
						}
					}

					polygon.erase(std::remove_if(polygon.begin(), polygon.end(), [](workingEdge& e) { return !e.valid_; }), polygon.end());

					std::for_each(polygon.begin(), polygon.end(), [&workingTriangles, &v](workingEdge& e) { workingTriangles.push_back(workingTriangle(e.i_, e.j_, &v)); });
				
				});

			// Remove the original super triangle...
			workingTriangles.erase(std::remove_if(workingTriangles.begin(), workingTriangles.end(), 
				[&p1Index, &p2Index, &p3Index](workingTriangle& t) 
				{ 
					return (t.containsVertex(p1Index) ||
						t.containsVertex(p2Index) ||
						t.containsVertex(p3Index));
				}), workingTriangles.end());

			// Create our edges (required for constraint functionality) and triangles...
			for (auto tItr = workingTriangles.begin(); tItr != workingTriangles.end(); ++tItr)
			{
				triangles_.push_back(triangle{ tItr->a_->index_, tItr->b_->index_, tItr->c_->index_ });
				edges_.push_back(edge{ tItr->a_->index_, tItr->b_->index_ });
				edges_.push_back(edge{ tItr->b_->index_, tItr->c_->index_ });
				edges_.push_back(edge{ tItr->c_->index_, tItr->a_->index_ });
			}
			
			// Remove our duplicate edges...
			removeDuplicateEdges();
		}

		void removeDuplicateEdges()
		{
			std::for_each(edges_.begin(), edges_.end(),
				[](edge& e)
				{
					if (e.j_ < e.i_)
					{
						unsigned int j = e.j_;
						e.j_ = e.i_;
						e.i_ = j;
					}
				});
			edges_.sort([](edge& e, edge& ee)
				{
					return e.i_ == ee.i_ ? e.j_ < ee.j_ : e.i_ < ee.i_;
				});
			edges_.unique([](edge& e, edge& ee)
				{
					return (e.i_ == ee.i_ && e.j_ == ee.j_) || (e.i_ == ee.j_ && e.j_ == ee.i_);
				});
		}

		std::list<triangle> triangles_;
		std::list<edge> edges_;
		std::vector<vector2> vertices_;
	};

}	// namespace delaunay








//
//
//namespace Delaunay
//{
//	const static double pi = 3.14159265358979323846;
//
//	static bool almost_equal(const double& x, const double& y, int ulp = 2)
//	{
//		// the machine epsilon has to be scaled to the magnitude of the values used
//		// and multiplied by the desired precision in ULPs (units in the last place)
//		return std::abs(x - y) <= std::numeric_limits<double>::epsilon() * std::abs(x + y) * static_cast<double>(ulp)
//			// unless the result is subnormal
//			|| std::abs(x - y) < (std::numeric_limits<double>::min)();
//	}
//
//	static double half(const double x)
//	{
//		return 0.5 * x;
//	}
//
//	struct Vector2
//	{
//		Vector2() = default;
//		Vector2(const Vector2 &v) = default;
//		Vector2(Vector2&&) = default;
//		Vector2(const double vx, const double vy)
//			: x(vx), y(vy)
//		{
//		}
//
//		double dist2(const Vector2 &v) const
//		{
//			const double dx = x - v.x;
//			const double dy = y - v.y;
//			return dx * dx + dy * dy;
//		}
//
//		double dist(const Vector2 &v) const
//		{
//			return hypot(x - v.x, y - v.y);
//		}
//
//		double norm2() const
//		{
//			return x * x + y * y;
//		}
//
//		Vector2 &operator=(const Vector2& v) = default;
//
//		Vector2 &operator=(Vector2&&) = default;
//
//		bool operator ==(const Vector2 &v) const
//		{
//			return (x == v.x) && (y == v.y);
//		}
//
//		double x;
//		double y;
//	};
//
//	static bool almost_equal(const Vector2 &v1, const Vector2 &v2, int ulp = 2)
//	{
//		return almost_equal(v1.x, v2.x, ulp) && almost_equal(v1.y, v2.y, ulp);
//	}
//
//	struct Edge
//	{
//		using VertexType = Vector2;
//
//		Edge() = default;
//		Edge(const Edge&) = default;
//		Edge(Edge&&) = default;
//
//		Edge &operator=(const Edge&) = default;
//		Edge &operator=(Edge&&) = default;
//
//		Edge(const VertexType &v1, const VertexType &v2)
//			: v(&v1), w(&v2)
//		{
//		}
//
//		Edge(VertexType* v1, VertexType* v2)
//			: v(v1), w(v2)
//		{
//		}
//
//		bool operator==(const Edge &e) const
//		{
//			return (*(this->v) == *e.v && *(this->w) == *e.w) ||
//				(*(this->v) == *e.w && *(this->w) == *e.v);
//		}
//		
//		const VertexType *v;
//		const VertexType *w;
//		bool isBad = false;
//	};
//
//	static bool almost_equal(const Edge &e1, const Edge &e2)
//	{
//		return	(almost_equal(*e1.v, *e2.v) && almost_equal(*e1.w, *e2.w)) ||
//			(almost_equal(*e1.v, *e2.w) && almost_equal(*e1.w, *e2.v));
//	}
//
//	struct Triangle
//	{
//		using EdgeType = Edge;
//		using VertexType = Vector2;
//
//		Triangle() = default;
//		Triangle(const Triangle&) = default;
//		Triangle(Triangle&&) = default;
//		Triangle(const VertexType &v1, const VertexType &v2, const VertexType &v3)
//			: a(&v1), b(&v2), c(&v3), isBad(false)
//		{}
//
//		bool containsVertex(const VertexType &v) const
//		{
//			return almost_equal(*a, v) || almost_equal(*b, v) || almost_equal(*c, v);
//		}
//
//		bool constainsEdge(const Edge& e) const
//		{
//			return containsVertex(*e.v) && containsVertex(*e.w);
//		}
//
//		bool circumCircleContains(const VertexType &v) const
//		{
//			const double ab = a->norm2();
//			const double cd = b->norm2();
//			const double ef = c->norm2();
//
//			const double ax = a->x;
//			const double ay = a->y;
//			const double bx = b->x;
//			const double by = b->y;
//			const double cx = c->x;
//			const double cy = c->y;
//
//			const double circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
//			const double circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));
//
//			const VertexType circum(half(circum_x), half(circum_y));
//			const double circum_radius = a->dist2(circum);
//			const double dist = v.dist2(circum);
//			return dist <= circum_radius;
//		}
//
//		Triangle &operator=(const Triangle&) = default;
//		Triangle &operator=(Triangle&&) = default;
//
//		bool operator ==(const Triangle &t) const
//		{
//			return	(*this->a == *t.a || *this->a == *t.b || *this->a == *t.c) &&
//				(*this->b == *t.a || *this->b == *t.b || *this->b == *t.c) &&
//				(*this->c == *t.a || *this->c == *t.b || *this->c == *t.c);
//		}
//
//		const VertexType *a;
//		const VertexType *b;
//		const VertexType *c;
//		bool isBad = false;
//	};
//
//	static bool almost_equal(const Triangle &t1, const Triangle &t2)
//	{
//		return	(almost_equal(*t1.a, *t2.a) || almost_equal(*t1.a, *t2.b) || almost_equal(*t1.a, *t2.c)) &&
//			(almost_equal(*t1.b, *t2.a) || almost_equal(*t1.b, *t2.b) || almost_equal(*t1.b, *t2.c)) &&
//			(almost_equal(*t1.c, *t2.a) || almost_equal(*t1.c, *t2.b) || almost_equal(*t1.c, *t2.c));
//	}
//
//	static double angleDifference(const Vector2& a, const Vector2& b)
//	{
//		double angle = atan2(b.y, b.x) - atan2(a.y, a.x);
//
//		if (angle > pi)
//			angle -= 2 * pi;
//		else if (angle <= -pi)
//			angle += 2 * pi;
//
//		return angle;
//	}
//
//	static std::optional<Vector2> segmentSegmentIntersection(const Vector2& a, const Vector2& b, const Vector2& i, const Vector2& j)
//	{
//		double s1_x, s1_y, s2_x, s2_y;
//		s1_x = b.x - a.x;     s1_y = b.y - a.y;
//		s2_x = j.x - i.x;     s2_y = j.y - i.y;
//
//		double sDenom = (-s2_x * s1_y + s1_x * s2_y);
//		double tDenom = (-s2_x * s1_y + s1_x * s2_y);
//
//		if (!sDenom || !tDenom)
//			return std::optional<Vector2>();
//
//		double s = (-s1_y * (a.x - i.x) + s1_x * (a.y - i.y)) / sDenom;
//		double t = (s2_x * (a.y - i.y) - s2_y * (a.x - i.x)) / tDenom;
//
//		if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
//			return std::optional<Vector2>(Vector2(a.x + (t * s1_x), a.y + (t * s1_y)));
//		
//		return std::optional<Vector2>();
//	}
//
//	class Tessellator
//	{
//	public:
//		using TriangleType = Triangle;
//		using EdgeType = Edge;
//		using VertexType = Vector2;
//
//		const std::list<TriangleType>& triangulate(const std::vector<VertexType> &vertices)
//		{
//			// Store the vertices locally
//			vertices_ = vertices;
//
//			// Determinate the super triangle
//			double minX = vertices_[0].x;
//			double minY = vertices_[0].y;
//			double maxX = minX;
//			double maxY = minY;
//
//			for (std::size_t i = 0; i < vertices_.size(); ++i)
//			{
//				if (vertices_[i].x < minX) minX = vertices_[i].x;
//				if (vertices_[i].y < minY) minY = vertices_[i].y;
//				if (vertices_[i].x > maxX) maxX = vertices_[i].x;
//				if (vertices_[i].y > maxY) maxY = vertices_[i].y;
//			}
//
//			const double dx = maxX - minX;
//			const double dy = maxY - minY;
//			const double deltaMax = (std::max)(dx, dy);
//			const double midx = half(minX + maxX);
//			const double midy = half(minY + maxY);
//
//			const VertexType p1(midx - 20 * deltaMax, midy - deltaMax);
//			const VertexType p2(midx, midy + 20 * deltaMax);
//			const VertexType p3(midx + 20 * deltaMax, midy - deltaMax);
//
//			// Create a list of triangles, and add the supertriangle in it
//			triangles_.push_back(TriangleType(p1, p2, p3));
//
//			for (auto p = begin(vertices_); p != end(vertices_); p++)
//			{
//				std::vector<EdgeType> polygon;
//
//				for (auto & t : triangles_)
//				{
//					if (t.circumCircleContains(*p))
//					{
//						t.isBad = true;
//						polygon.push_back(Edge{ *t.a, *t.b });
//						polygon.push_back(Edge{ *t.b, *t.c });
//						polygon.push_back(Edge{ *t.c, *t.a });
//					}
//				}
//
//				triangles_.erase(std::remove_if(begin(triangles_), end(triangles_), [](TriangleType &t) {
//					return t.isBad;
//				}), end(triangles_));
//
//				for (auto e1 = begin(polygon); e1 != end(polygon); ++e1)
//				{
//					for (auto e2 = e1 + 1; e2 != end(polygon); ++e2)
//					{
//						if (almost_equal(*e1, *e2))
//						{
//							e1->isBad = true;
//							e2->isBad = true;
//						}
//					}
//				}
//
//				polygon.erase(std::remove_if(begin(polygon), end(polygon), [](EdgeType &e) {
//					return e.isBad;
//				}), end(polygon));
//
//				for (const auto e : polygon)
//					triangles_.push_back(TriangleType(*e.v, *e.w, *p));
//			}
//
//			triangles_.erase(std::remove_if(begin(triangles_), end(triangles_), [p1, p2, p3](TriangleType &t) {
//				return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
//			}), end(triangles_));
//
//			for (const auto t : triangles_)
//			{
//				edges_.push_back(Edge{ *t.a, *t.b });
//				edges_.push_back(Edge{ *t.b, *t.c });
//				edges_.push_back(Edge{ *t.c, *t.a });
//			}
//
//			// remove duplicate edges...
//			for (auto eItr = edges_.begin(); eItr != edges_.end(); ++eItr)
//			{
//
//			}
//
//
//			unsigned int edgeCount = edges_.size();
//			auto eItr = edges_.begin();
//			for (unsigned int e = 0; e < edgeCount; ++e, ++eItr)
//				edges_.push_back(Edge{*eItr->w, *eItr->v});
//			
//			edges_.sort([](EdgeType& e, EdgeType& ee) { return e.v == ee.v ? e.w < ee.w : e.v < ee.v; });
//			edges_.unique([](EdgeType& e, EdgeType& ee) { return almost_equal(ee, e); });
//			
//			return triangles_;
//		}
//
//		const std::list<TriangleType>& getTriangles() const { return triangles_; }
//		const std::list<EdgeType>& getEdges() const { return edges_; }
//		const std::vector<VertexType>& getVertices() const { return vertices_; }
//
//		void addConstraint(const std::pair<unsigned int, unsigned int>& e)
//		{
//			// check if we already have this edge...
//			VertexType* e1 = &vertices_[e.first];
//			VertexType* e2 = &vertices_[e.second];
//			
//			if (std::find_if(edges_.begin(), edges_.end(), [&e1, &e2](const EdgeType& edge) { return ((e1 == edge.v && e2 == edge.w) || (e1 == edge.w && e2 == edge.v)); }) == edges_.end())
//			{
//				Vector2 constraintEdge(e2->x - e1->x, e2->y - e1->y);
//
//				// find out which edges are intersected, and for the two vertices of the intersected edge, find out which side of the constraint they are on...
//				std::list<VertexType> leftSide, rightSide;
//				for (std::list<EdgeType>::iterator ee = edges_.begin(); ee != edges_.end();)
//				{
//					if (*ee->v == *e1 || *ee->w == *e1 || *ee->v == *e2 || *ee->w == *e2)
//					{
//						++ee;
//						continue;
//					}
//						
//
//					double vAngle = angleDifference(Vector2(ee->v->x - e1->x, ee->v->y - e1->y), constraintEdge);
//					double wAngle = angleDifference(Vector2(ee->w->x - e2->x, ee->w->y - e2->y), constraintEdge);
//
//					if ((vAngle < 0 && wAngle > 0) || (vAngle > 0 && wAngle < 0))
//					{
//						// do the addtional check to see if this constrain edge is intersected...
//						if (segmentSegmentIntersection(*e1, *e2, *ee->v, *ee->w))
//						{
//							if (vAngle > 0)
//								leftSide.push_back(*ee->v);
//							else
//								rightSide.push_back(*ee->v);
//
//							if (wAngle > 0)
//								leftSide.push_back(*ee->w);
//							else
//								rightSide.push_back(*ee->w);
//
//							// and any triangles that use this edge...
//							triangles_.erase(std::remove_if(begin(triangles_), end(triangles_), [&ee](Triangle& tri) { return tri.constainsEdge(*ee); }), end(triangles_));
//
//							// remove the edge from our list...
//							ee = edges_.erase(ee);
//
//							// next edge...
//							continue;
//						}
//					}
//
//					++ee;
//				}
//
//				// add the constaint edge to both sides and remove the duplicates of both sides.
//				leftSide.push_back(*e1);
//				leftSide.push_back(*e2);
//				rightSide.push_back(*e1);
//				rightSide.push_back(*e2);
//
//				// tessellate each subset...
//				Tessellator leftTessellation, rightTessellation;
//
//				// remove duplicates and tessellation the left side...
//				leftSide.sort([](const VertexType& a, const VertexType& b) { if (a.x == b.x) { return a.y < b.y; } else { return a.x < b.x; } });
//				leftSide.erase(std::unique(leftSide.begin(), leftSide.end()), leftSide.end());
//				leftTessellation.triangulate(std::vector<VertexType>(leftSide.begin(), leftSide.end()));
//
//				// remove duplicates and tessellation the right side...
//				rightSide.sort([](const VertexType& a, const VertexType& b) { if (a.x == b.x) { return a.y < b.y; } else { return a.x < b.x; } });
//				rightSide.erase(std::unique(rightSide.begin(), rightSide.end()), rightSide.end());
//				rightTessellation.triangulate(std::vector<VertexType>(rightSide.begin(), rightSide.end()));
//
//				// add these two subsets to the final results...
//				append(leftTessellation);
//				append(rightTessellation);
//			}
//		}
//
//	private:
//		VertexType* addVertex(const VertexType& v) { vertices_.push_back(v); return &vertices_.back(); }
//
//		// destructive copy...
//		void append(Tessellator& tessellation)
//		{
//			// create new container...
//			std::vector<VertexType> newVertices;
//			
//			// copy the vertices across...
//			newVertices.insert(newVertices.end(), vertices_.begin(), vertices_.end());
//			newVertices.insert(newVertices.end(), tessellation.vertices_.begin(), tessellation.vertices_.end());
//
//			// iterate through the original vertices of this tessellation and set the x value to the index of it.
//			for (unsigned int v = 0; v < vertices_.size(); ++v)
//				vertices_[v].x = static_cast<double>(v);
//
//			// iterate through the original vertices of the tessellation to append and set the x value to the index of it + offset of original this tessellation vertices.
//			unsigned int offset = vertices_.size();
//			for (unsigned int v = 0; v < tessellation.vertices_.size(); ++v)
//				tessellation.vertices_[v].x = static_cast<double>(v + offset);
//
//			// set all the correct edge triangle pointers. new offset of new container is supplied by 
//
//			// now copy the edges and triangles across to the new container, set the pointer to the new index supplied by x coord.
//			triangles_.insert(triangles_.end(), tessellation.triangles_.begin(), tessellation.triangles_.end());
//			std::for_each(triangles_.begin(), triangles_.end(), 
//				[&newVertices](TriangleType& t) 
//			{
//				t.a = &newVertices[static_cast<unsigned int>(t.a->x)];
//				t.b = &newVertices[static_cast<unsigned int>(t.b->x)];
//				t.c = &newVertices[static_cast<unsigned int>(t.c->x)];
//			});
//			
//			// do the same for the edges...
//			edges_.insert(edges_.end(), tessellation.edges_.begin(), tessellation.edges_.end());
//			std::for_each(edges_.begin(), edges_.end(),
//				[&newVertices](EdgeType& e)
//			{
//				e.v = &newVertices[static_cast<unsigned int>(e.v->x)];
//				e.w = &newVertices[static_cast<unsigned int>(e.w->x)];
//			});
//
//			// swap vertices...
//			vertices_.swap(newVertices);
//		}
//
//		std::list<TriangleType> triangles_;
//		std::list<EdgeType> edges_;
//		std::vector<VertexType> vertices_;
//	};
//
//
//	std::string csv(const std::vector<Vector2>& points)
//	{
//		std::ostringstream oss;
//		std::for_each(points.begin(), points.end(), 
//			[&oss](const Vector2& v) 
//		{ 
//			oss << v.x << ", " << v.y << '\n'; 
//		});
//		return oss.str();
//	}
//
//
//	void validate(std::vector<Vector2>& vertices)
//	{
//		std::sort(vertices.begin(), vertices.end(), 
//			[](const Vector2& a, const Vector2& b) 
//		{
//			if (a.x == b.x)
//				return a.y < b.y;
//			return a.x < b.x;
//		});
//
//		vertices.erase(std::unique(vertices.begin(), vertices.end()));
//	}
//}

#endif // DELAUNAY_TRIANGULATION_HPP