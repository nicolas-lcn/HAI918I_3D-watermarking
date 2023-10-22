#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "Vec3.h"

class Triangle;


class Mesh 
{
public:
    std::vector< Vec3 > vertices; //array of mesh vertices positions
    std::vector< Vec3 > normals; //array of vertices normals useful for the display
    std::vector< Triangle > triangles; //array of mesh triangles
    std::vector< Vec3 > triangle_normals; //triangle normals to display face normals
    std::vector<Vec3> boardingbox;

public:

	Mesh();
	~Mesh();


public:
    //Compute face normals for the display
    void computeTrianglesNormals();

    //Compute vertices normals as the average of its incident faces normals
    void computeVerticesNormals(int type);

    void computeNormals(int type);

    //return boarding box of the mesh
    void boardingBox();

    Vec3 barycentre();

    // returns min and max radial coordinate of vertices
    std::vector<float> radialMinMax();

    std::vector<float> normalize(std::vector<float> i_values, float min, float max);
    std::vector<float> unnormalize(std::vector<float> normalizedValues, float min, float max);

    

};


#endif