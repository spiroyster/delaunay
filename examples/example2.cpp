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

#include "..\include\delaunay.hpp"

#include <fstream>
#include <iostream>
#include <exception>

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

std::string writeOBJ(const delaunay::tessellator& tessellator)
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
			oss << "f " << t.a_ + 1 << " " << t.b_ + 1 << " " << t.c_ + 1 << '\n';
		});

	return oss.str();
}

int main(int argc, char** argv)
{
	try
	{
		std::string verticesFilename;
		std::string edgesFilename;

		// read the arguments...
		if (argc == 1)
			throw std::exception("No input files.");
		if (argc >= 2)
			verticesFilename = std::string(argv[1]);
		if (argc == 3)
			edgesFilename = std::string(argv[2]);
		if (argc > 3)
			throw std::exception("Too many arguments.");

		// read the csv vertices...
		std::vector<delaunay::vector2> inputVertices = readCSV(verticesFilename);
		std::vector<delaunay::vector2> inputEdges = readCSV(edgesFilename);

		// perform tessellation...
		delaunay::tessellator dt(inputVertices);

		// add constraints...
		std::for_each(inputEdges.begin(), inputEdges.end(), [&dt](const delaunay::vector2& edge) { dt.addConstraint(delaunay::edge{ static_cast<unsigned int>(edge.x_), static_cast<unsigned int>(edge.y_) }); });

		// output obj...
		std::cout << writeOBJ(dt);

		return 0;
	}
	catch (const std::exception& e)
	{
		std::cout << "-= ERROR =- " << e.what() << "\n\n";

		// output help message...
		std::cout << "(Bowyer-Watson) Constrained Delaunay Tessellator - command line program.\n";
		std::cout << "https://github.com/spiroyster/delaunay/examples/example2.cpp\n";
		std::cout << "\n"; 
		std::cout << "Usage:   example2 {vertices.csv} {edges.csv} > output.obj\n";
		std::cout << "\n";
		std::cout << "Input Vertices: 2D csv values (double, double).\n";
		std::cout << "Input Edges: 2D csv values (int, int).\n";
		std::cout << "Output: stdout (wavefront) obj format.\n";
	}
	
	return 1;
}
