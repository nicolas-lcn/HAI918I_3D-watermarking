#include "VoxelGrid.h"

#include "Mesh.h"


VoxelGrid::VoxelGrid(Vec3 _minBound, Vec3 _maxBound, int _resolution):
	minBound(_minBound), 
	maxBound(_maxBound), 
	grid(std::map<int, Voxel>()), 
	nx(_resolution), 
	ny(_resolution), 
	nz(_resolution)
{
	dx = (maxBound[0] - minBound[0])/nx;
	dy = (maxBound[1] - minBound[1])/ny;
	dz = (maxBound[2] - minBound[2])/nz;
}

VoxelGrid::~VoxelGrid(){}

size_t VoxelGrid::size() { return grid.size();}

void VoxelGrid::computeVoxel(Vec3 vertice, Vec3 normal)
{
	int id = getGridIndex(toGridCoordinates(vertice));
	if(grid.count(id) != 0)
	{
		grid[id].pos += vertice;
		grid[id].normal += normal;
		grid[id].nbVertices ++;
	}
	else
	{
		Voxel voxel;
		voxel.pos = vertice;
		voxel.normal = normal;
		voxel.nbVertices = 1;
		grid[id] = voxel;
	}
	
}

void VoxelGrid::normalizeAll()
{
	 for (auto& [id, voxel] : grid)
	 {
	 	voxel.pos /= voxel.nbVertices;
		voxel.normal.normalize();
	 }
}

int VoxelGrid::getGridIndex(Vec3 gridPos)
{
	return gridPos[0] + gridPos[1] * nx + gridPos[2] * nx * ny;
}

Vec3 VoxelGrid::toGridCoordinates(Vec3 pos)
{
	return Vec3((int)((pos[0]-minBound[0])/dx), (int)((pos[1]-minBound[1])/dy), (int)((pos[2]-minBound[2])/dz));
}

Vec3 VoxelGrid::toWorldCoordinates(Vec3 gridPos)
{
	return Vec3(gridPos[0]*dx+minBound[0], gridPos[1]*dy+minBound[1], gridPos[2]*dz+minBound[2]);
}

Voxel VoxelGrid::getVoxelAt(int index)
{
	return grid[index];
}