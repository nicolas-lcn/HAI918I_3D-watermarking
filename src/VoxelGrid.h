#ifndef VOXELGRID_H
#define VOXELGRID_H

#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "Vec3.h"

class Mesh;

//Representant
struct Voxel
{
	Vec3 pos;
	Vec3 normal;
	int nbVertices; //nb vertices this struct is representing
};

class VoxelGrid
{
private:

	//bounds in world coordinates;
	Vec3 minBound;
	Vec3 maxBound;
	float dx, dy, dz;
	float nx, ny, nz;

	std::map<int, Voxel> grid;

public:

	size_t size();

	VoxelGrid(Vec3 minBound, Vec3 maxBound, int resolution);
	~VoxelGrid();

	int getGridIndex(Vec3 gridPos);

	Vec3 toGridCoordinates(Vec3 pos);

	Vec3 toWorldCoordinates(Vec3 gridPos);

	void computeVoxel(Vec3 vertice, Vec3 normal);

	Voxel getVoxelAt(int index);

	void normalizeAll();

};


#endif