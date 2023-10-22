#ifndef COMPRESSOR_H
#define COMPRESSOR_H
#include <cmath>
#include <algorithm>
#include "Vec3.h"


class Mesh;

class Compressor
{
public:
	static void quantify(Mesh &mesh, int qp)
	{
		float range;
		std::vector<Vec3> bb = mesh.boardingbox;
		Vec3 bbmin = bb[0];
		range = std::max((bb[1] - bb[0]).length(), std::max((bb[2] - bb[0]).length(), (bb[3] - bb[0]).length()));
		for (int i = 0; i < mesh.vertices.size(); ++i)
		{
			int x,y,z;
			x = (int) mesh.vertices[i][0];		
			y = (int) mesh.vertices[i][1];		
			z = (int) mesh.vertices[i][2];		
			int x_p, y_p, z_p;
			x_p = (x - (int)bbmin[0])*pow(2, qp)/range;
			y_p = (y - (int)bbmin[1])*pow(2, qp)/range;
			z_p = (z - (int)bbmin[2])*pow(2, qp)/range;
			mesh.vertices[i] = Vec3((float)x_p, (float)y_p, (float)z_p);

		}	
	}

	static void dequantify(Mesh &mesh, int qp)
	{
		float range;
		std::vector<Vec3> bb = mesh.boardingbox;
		Vec3 bbmin = bb[0];
		range = std::max((bb[1] - bb[0]).length(), std::max((bb[2] - bb[0]).length(), (bb[3] - bb[0]).length()));
		for (int i = 0; i < mesh.vertices.size(); ++i)
		{
			Vec3 c_p = mesh.vertices[i];
			Vec3 c = range/pow(2, qp) * c_p + bbmin;
			mesh.vertices[i] = c;
		}
	}

	static std::vector<Vec3> quant_dequant(Mesh &mesh, int qp)
	{
		float range;
		std::vector<Vec3> bb = mesh.boardingbox;
		Vec3 bbmin = bb[0];
		std::vector<Vec3> result = std::vector<Vec3>(mesh.vertices.size());
		range = std::max((bb[1] - bb[0]).length(), std::max((bb[2] - bb[0]).length(), (bb[3] - bb[0]).length()));
		for (int i = 0; i < mesh.vertices.size(); ++i)
		{
			float x,y,z;
			x = mesh.vertices[i][0];		
			y = mesh.vertices[i][1];		
			z = mesh.vertices[i][2];		
			int x_p, y_p, z_p;
			x_p = (int)((x - bbmin[0])*pow(2, qp)/range);
			y_p = (int)((y - bbmin[1])*pow(2, qp)/range);
			z_p = (int)((z - bbmin[2])*pow(2, qp)/range);
			float x_p2, y_p2, z_p2;
			x_p2 = (float)(x_p*range/pow(2,qp)) + bbmin[0];
			y_p2 = (float)(y_p*range/pow(2,qp)) + bbmin[1];
			z_p2 = (float)(z_p*range/pow(2,qp)) + bbmin[2];

			result[i] = Vec3(x_p2, y_p2, z_p2);

		}
		return result;
		
	}

	static double RMSE(std::vector<Vec3> v0, std::vector<Vec3> v1)
	{
		// sqrt (1/N * sum 1 à N de y0 - y1 ²)
		float sum = 0;
		for (int i = 0; i < v0.size(); ++i)
		{
			sum += pow((v0[i][0] - v1[i][0]), 2) + pow((v0[i][1] - v1[i][1]), 2) + pow((v0[i][2] - v1[i][2]), 2);
		}
		sum/=(float)v0.size();
		float rmse = sqrt(sum);
		return rmse;
		
	}	
};


#endif 