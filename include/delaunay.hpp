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

	part of https://github.com/spiroyster/delaunay
*/

#ifndef DELAUNAY_TESSELLATION_HPP
#define DELAUNAY_TESSELLATION_HPP

#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include <assert.h>
#include <string>
#include <sstream>

namespace delaunay
{
	struct vector2 { vector2(double x, double y) : x_(x), y_(y) {} double x_, y_; };
	struct edge { unsigned int i_, j_; };
	struct triangle { unsigned int a_, b_, c_; };

	const static double pi = 3.14159265358979323846;

	namespace convenience
	{
		std::string csv(const std::vector<vector2>& points)
		{
			std::ostringstream oss;
			std::for_each(points.begin(), points.end(),[&oss](const vector2& v) { oss << v.x_ << ", " << v.y_ << '\n'; });
			return oss.str();
		}

		void validate(std::vector<vector2>& vertices)
		{
			std::sort(vertices.begin(), vertices.end(), [](const vector2& a, const vector2& b) { return a.x_ == b.x_ ? a.y_ < b.y_ : a.x_ < b.x_; });
			vertices.erase(std::unique(vertices.begin(), vertices.end(), [](const vector2& a, const vector2& b) { return a.x_ == b.x_ && a.y_ == b.y_; }), vertices.end());
		}
	}

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

		static bool orEquals(unsigned int a, unsigned int b, unsigned int c, unsigned int equals)
		{
			return (a == equals || b == equals || c == equals);
		}

		struct workingTriangle
		{
			workingTriangle(const workingVertex* a, const workingVertex* b, const workingVertex* c)
				: a_(a), b_(b), c_(c), valid_(true) {}

			bool containsVertex(unsigned int index)
			{
				return orEquals(a_->index_, b_->index_, c_->index_, index);
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

		double angleDifference(const vector2& a, const vector2& b)
		{
			double angle = atan2(b.y_, b.x_) - atan2(a.y_, a.x_);
			
			if (angle > pi)
				angle -= 2 * pi;
			else if (angle <= -pi)
				angle += 2 * pi;
			
			return angle;
		}

		bool segmentSegmentIntersection(const vector2& a, const vector2& b, const vector2& c, const vector2& d)
		{
			double s1_x, s1_y, s2_x, s2_y;
			s1_x = b.x_ - a.x_;     s1_y = b.y_ - a.y_;
			s2_x = d.x_ - c.x_;     s2_y = d.y_ - c.y_;
			
			double sDenom = (-s2_x * s1_y + s1_x * s2_y);
			double tDenom = (-s2_x * s1_y + s1_x * s2_y);
			
			if (!sDenom || !tDenom)
				return false;
			
			double s = (-s1_y * (a.x_ - c.x_) + s1_x * (a.y_ - c.y_)) / sDenom;
			double t = (s2_x * (a.y_ - c.y_) - s2_y * (a.x_ - c.x_)) / tDenom;
			
			return (s >= 0 && s <= 1 && t >= 0 && t <= 1);
		}
		

	public:
		tessellator(const std::vector<vector2>& points)
			: vertices_(points)
		{
			#ifdef NDEBUG
			std::vector<vector2> verticesCopy = points;
			convenience::validate(verticesCopy);
			assert(verticesCopy.size() == points.size());
			#endif

			std::vector<unsigned int> allVertexIndexes(points.size(), 0);
			for (unsigned int i = 0; i < points.size(); ++i)
				allVertexIndexes[i] = i;
			
			tessellate(allVertexIndexes);

			// Remove our duplicate edges...
			removeDuplicateEdges();
		}

		void addConstraint(const edge& constraint)
		{
			// assert the requested constraint indexes are valid...
			assert(constraint.i_ < vertices_.size() && constraint.j_ < vertices_.size());

			// Check if this edge is already present, if so do nothing...
			if (std::find_if(edges_.begin(), edges_.end(), [&constraint](const edge& ee) { return (constraint.i_ == ee.i_ && constraint.j_ == ee.j_) || (constraint.i_ == ee.j_ && constraint.j_ == ee.i_); }) != edges_.end())
				return;

			std::vector<unsigned int> leftSide, rightSide;

			// find out which edges this constraint straddles...
			for (auto e = edges_.begin(); e != edges_.end();)
			{
				// Don't check edges which share a constraint vertex...
				if (constraint.i_ != e->i_ && constraint.i_ != e->j_ && constraint.j_ != e->i_ && constraint.j_ != e->j_)
				{
					const vector2& ci = vertices_[constraint.i_];
					const vector2& cj = vertices_[constraint.j_];
					const vector2& ei = vertices_[e->i_];
					const vector2& ej = vertices_[e->j_];
					
					double iAngle = angleDifference(vector2(ei.x_ - ci.x_, ei.y_ - ci.y_), vector2(cj.x_ - ci.x_, cj.y_ - ci.y_));
					double jAngle = angleDifference(vector2(ej.x_ - ci.x_, ej.y_ - ci.y_), vector2(cj.x_ - ci.x_, cj.y_ - ci.y_));

					if ((iAngle < 0 && jAngle > 0) || (iAngle > 0 && jAngle < 0))
					{
						// constraint straddles an edge, so check for edge intersection...
						if (segmentSegmentIntersection(ei, ej, ci, cj))
						{
							leftSide.push_back(iAngle > 0 ? e->i_ : e->j_);
							rightSide.push_back(iAngle < 0 ? e->i_ : e->j_);
							
							// Remove any triangles using the edge, and the edge itself...
							triangles_.erase(std::remove_if(triangles_.begin(), triangles_.end(), [&e](triangle& t) { return orEquals(t.a_, t.b_, t.c_, e->i_) && orEquals(t.a_, t.b_, t.c_, e->j_); }), triangles_.end());
							e = edges_.erase(e);

							continue;
						}
					}
				}

				++e;
			}

			// add the constraint edge to both sides, and tessellate subsets...
			leftSide.push_back(constraint.i_);
			leftSide.push_back(constraint.j_);
			rightSide.push_back(constraint.i_);
			rightSide.push_back(constraint.j_);

			// remove duplicate indexes...
			std::sort(leftSide.begin(), leftSide.end());
			leftSide.erase(std::unique(leftSide.begin(), leftSide.end()), leftSide.end());
			std::sort(rightSide.begin(), rightSide.end());
			rightSide.erase(std::unique(rightSide.begin(), rightSide.end()), rightSide.end());

			// retessellate both the right and left side
			tessellate(leftSide);
			tessellate(rightSide);

			// Remove our duplicate edges...
			removeDuplicateEdges();
		}

		const std::list<triangle>& getTriangles() const { return triangles_; }
		const std::vector<vector2>& getVertices() const { return vertices_; }
		const std::list<edge>& getEdges() const { return edges_; }

	private:
		void tessellate(const std::vector<unsigned int>& vertexIndexes)
		{
			#ifdef NDEBUG
			assert(!vertexIndexes.empty());
			assert(!vertices_.empty());
			for (unsigned int vi = 0; vi < vertexIndexes.size(); ++vi)
				assert(vertexIndexes[vi] < vertices_.size());
			std::vector<unsigned int> vertexIndexesCopy = vertexIndexes;
			std::sort(vertexIndexesCopy.begin(), vertexIndexesCopy.end(), [](const unsigned int& a, const unsigned int& b) { return a < b; });
			assert(std::unique(vertexIndexesCopy.begin(), vertexIndexesCopy.end()) == vertexIndexesCopy.end())
			#endif

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
			unsigned int p1Index = vertices_.size(), p2Index = vertices_.size() + 1, p3Index = vertices_.size() + 2;

			vector2 p1(midx - 20.0 * deltaMax, midy - deltaMax);
			vector2 p2(midx, midy + 20.0 * deltaMax);
			vector2 p3(midx + 20.0 * deltaMax, midy - deltaMax);

			workingVertices[vertexIndexes.size()] = workingVertex(&p1, p1Index);
			workingVertices[vertexIndexes.size()+1] = workingVertex(&p2, p2Index);
			workingVertices[vertexIndexes.size()+2] = workingVertex(&p3, p3Index);

			// Create a list of triangles, and add the supertriangle in it
			std::list<workingTriangle> workingTriangles;
			workingTriangles.push_back(workingTriangle(&workingVertices[vertexIndexes.size()], &workingVertices[vertexIndexes.size()+1], &workingVertices[vertexIndexes.size()+2]));

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

#endif // DELAUNAY_TESSELLATION_HPP