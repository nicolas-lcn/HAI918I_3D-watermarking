#include "Mesh.h"
#include "Triangle.h"
#include "Vec3.h"

Mesh::Mesh(){}
Mesh::~Mesh(){}


//Compute face normals for the display
void Mesh::computeTrianglesNormals(){
    triangle_normals.clear();
    for(size_t i = 0; i < triangles.size(); i++)
    {
        Vec3 e_10 = vertices[triangles[i][1]] - vertices[triangles[i][0]];
        Vec3 e_20 = vertices[triangles[i][2]] - vertices[triangles[i][0]];
        Vec3 normal = Vec3::cross(e_10,e_20);
        normal.normalize();
        triangle_normals.push_back(normal);
    }
}

void Mesh::computeVerticesNormals(int weight_type){  
    normals.clear();
    
    if(weight_type == 0)
    {
        for(size_t i = 0; i < vertices.size(); i++)
        {
            normals.push_back(Vec3(0,0,0));
        }
        
        for(size_t i = 0; i<triangles.size(); i++)
        {
            for(int j = 0 ; j<3; j++)
            {
                normals[triangles[i][j]] += triangle_normals[i];
            }
        }         

    }
    else if (weight_type == 1)
    {

        std::vector<float> areas;
        areas.resize(vertices.size());
        for(size_t i = 0; i < triangles.size(); i++)
        {
            Vec3 e_10 = vertices[triangles[i][1]] - vertices[triangles[i][0]];
            Vec3 e_20 = vertices[triangles[i][2]] - vertices[triangles[i][0]];
            float area = (Vec3::cross(e_10,e_20).length())/2.0f;
            areas.push_back(area);
        }

        for(size_t i = 0; i < vertices.size(); i++)
        {
            normals.push_back(Vec3(0,0,0));
        }
        
        for(size_t i = 0; i<triangles.size(); i++)
        {
            for(int j = 0 ; j<3; j++)
            {
                normals[triangles[i][j]] += areas[i] * triangle_normals[i];
            }
        } 
    }

    else
    {
        std::vector<float> angles;
        

        for(size_t i = 0; i < triangles.size(); i++)
        {
            Vec3 e_10 = vertices[triangles[i][1]] - vertices[triangles[i][0]];
            Vec3 e_20 = vertices[triangles[i][2]] - vertices[triangles[i][0]];
            float angle = acos(abs(Vec3::dot(e_10,e_20))/(e_10.length()*e_20.length()));
            angles.push_back(angle);
        }
        for(size_t i = 0; i < vertices.size(); i++)
        {
            normals.push_back(Vec3(0,0,0));
        }
        
        for(size_t i = 0; i<triangles.size(); i++)
        {
            for(int j = 0 ; j<3; j++)
            {
                normals[triangles[i][j]] += angles[i] * triangle_normals[i];
            }
        } 
    }
   
    for(Vec3 & normal: normals)
    {
        normal.normalize();
    }

}

void Mesh::computeNormals(int weight_type){
    computeTrianglesNormals();
    computeVerticesNormals(weight_type);
}

void Mesh::boardingBox()
{
	Vec3 bbmin, bbmax;
	bbmin = vertices[0];
	bbmax = vertices[0];

	for(size_t i = 0; i<vertices.size(); i++)
	{
		for(int j = 0; j<3; j++)
		{
			if(bbmin[j] > vertices[i][j]) bbmin[j] = vertices[i][j];
			if(bbmax[j] < vertices[i][j]) bbmax[j] = vertices[i][j];
		}
	}
	bbmin -= Vec3(0.01, 0.01, 0.01);
	bbmax += Vec3(0.01, 0.01, 0.01);

	boardingbox.push_back(bbmin);
	boardingbox.push_back(Vec3(bbmax[0], bbmin[1], bbmin[2]));
	boardingbox.push_back(Vec3(bbmin[0], bbmax[1], bbmin[2]));
	boardingbox.push_back(Vec3(bbmin[0], bbmin[1], bbmax[2]));
	boardingbox.push_back(Vec3(bbmax[0], bbmax[1], bbmin[2]));
	boardingbox.push_back(Vec3(bbmax[0], bbmin[1], bbmax[2]));
	boardingbox.push_back(Vec3(bbmin[0], bbmax[1], bbmax[2]));
	boardingbox.push_back(bbmax);

}

Vec3 Mesh::barycentre()
{
    Vec3 barycentre = Vec3(0,0,0);
    for (int i = 0; i < vertices.size(); ++i)
    {
        barycentre += vertices[i];
    }
    return barycentre/vertices.size();
}

std::vector<float> Mesh::radialMinMax()
{
    float min = FLT_MAX;
    float max = -FLT_MAX;
    for (int i = 0; i < vertices.size(); ++i)
    {
        Vec3 spherical = Vec3::EuclideanCoordinatesToSpherical(vertices[i]);
        if(spherical[2] <= min) min = spherical[2];
        if(spherical[2] >= max) max = spherical[2];
    }
    return {min, max};
}


