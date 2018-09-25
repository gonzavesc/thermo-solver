#include <iostream>
#include <fstream>
#include <vector>
#include <string> 
#include <cmath>
#include <new>

#ifndef INCLUDE_MATERIAL
	#include "material.hpp"
	#define INCLUDE_MATERIAL
#endif

#ifndef INCLUDE_TEST
	#include "test.hpp"
	#define INCLUDE_TEST
#endif

std::vector<double> readfiledat()
{
    std::ifstream file;
    std::vector<double> v;

    file.open("config.dat");    
    std::string name;
    while (file >> name) 
    {
        double a;

        file >> a;
        v.push_back(a);
    }
    return v;
}

void setvalues(const std::vector<double>& C, std::vector<Material>& M)
{
    M[0].rho = C[0]; M[0].cp = C[1]; M[0].k = C[2];
    M[1].rho = C[3]; M[1].cp = C[4]; M[1].k = C[5];
    M[2].rho = C[6]; M[2].cp = C[7]; M[2].k = C[8];
    M[3].rho = C[9]; M[3].cp = C[10]; M[3].k = C[11];
}

void discret_points(const std::vector<double>& delta, const std::vector<double>& C, std::vector<std::vector<int>>& points)
{
    points[0][0] = C[12] / delta[0]; points [1][0] = C[13] / delta [1];
    points[0][1] = C[14] / delta[0]; points [1][1] = C[15] / delta [1];
    points[0][2] = C[16] / delta[0]; points [1][2] = C[17] / delta [1];
}

void setcoeff(const std::vector<double>& diff, std::vector<Material>& M)
{
    int i;
    for (i = 0; i <= 3; i++)
    {
        M[i].ae = M[i].k * diff[1] / diff[0]; M[i].aw = M[i].ae;
        M[i].an = M[i].k * diff[0] / diff[1]; M[i].as = M[i].an;
        M[i].a0 = M[i].rho * M[i].cp * diff[0] * diff[1] / diff[2];
        M[i].ap = M[i].ae + M[i].aw + M[i].an + M[i].as + M[i].a0;
    }
}