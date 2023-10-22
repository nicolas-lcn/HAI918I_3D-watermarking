#ifndef SIMPLIFIER_H
#define SIMPLIFIER_H

class Mesh;

class Simplifier
{
public:
	static void simplify (unsigned int resolution, const Mesh &mesh, Mesh &result);
};

#endif