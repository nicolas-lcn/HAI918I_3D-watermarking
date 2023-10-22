#include "Simplifier.h"
#include "VoxelGrid.h"
#include "Mesh.h"
#include "Triangle.h"

void Simplifier::simplify(unsigned int resolution, const Mesh &mesh, Mesh &result)
{
	std::vector<Vec3> bb = mesh.boardingbox;
	VoxelGrid grid(bb.at(0), bb.at(bb.size()-1), resolution);

	for(size_t i = 0; i<mesh.vertices.size(); i++)
		grid.computeVoxel(mesh.vertices[i], mesh.normals[i]);
	grid.normalizeAll();
	size_t nbTriangles = mesh.triangles.size();
	//std::vector<size_t> remove;
	unsigned int count = 0;
	for(size_t i = 0; i<nbTriangles; i++)
	{
		int r0 = grid.getGridIndex(grid.toGridCoordinates(mesh.vertices[mesh.triangles[i][0]]));
		int r1 = grid.getGridIndex(grid.toGridCoordinates(mesh.vertices[mesh.triangles[i][1]]));
		int r2 = grid.getGridIndex(grid.toGridCoordinates(mesh.vertices[mesh.triangles[i][2]]));

		if(r0 == r1 || r0 == r2 || r1 == r2) continue;
		else
		{
			Voxel representant = grid.getVoxelAt(r0);
			result.vertices.push_back(representant.pos);
			unsigned int c0 = result.vertices.size()-1;
			result.normals.push_back(representant.normal);
			representant = grid.getVoxelAt(r1);
			result.vertices.push_back(representant.pos);
			unsigned int c1 = result.vertices.size()-1;
			result.normals.push_back(representant.normal);
			representant = grid.getVoxelAt(r2);
			result.vertices.push_back(representant.pos);
			unsigned int c2 = result.vertices.size()-1;
			result.normals.push_back(representant.normal);
			result.triangles.push_back(Triangle(c0, c1, c2));
			
		}
	}
}