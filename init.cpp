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

#ifndef INCLUDE_INIT
	#include "init.hpp"
	#define INCLUDE_INIT
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

void setcontact(const std::vector<double>& diff, Material& contact, const Material& M1, const Material& M2, const int dir, const int side)
{
    double keq;
    contact.ae = M1.ae; contact.aw = M1.aw; contact.an = M1.an; contact.as = M1.an;    
    contact.k = M1.k; contact.rho = M1.rho; contact.cp = M1.cp; contact.a0 = M1.a0;
    keq = 2 * M1.k * M2.k / (M1.k + M2.k);
    if (dir == 0)
    {
        if (side == 0)
        {
            contact.ae = keq * diff[1] / diff[0];
        }
        if (side == 1)
        {
            contact.aw = keq * diff[1] / diff[0];
        }
    }
    if (dir == 1)
    {
        if (side == 0)
        {
            contact.an = keq * diff[0] / diff[1];
        }
        if (side == 1)
        {
            contact.as = keq * diff[0] / diff[1];
        }
    }
    contact.ap = contact.ae + contact.aw + contact.an + contact.as + contact.a0;
}

void initializeT(std::vector<std::vector<double>>& Tmap, const int p)
{
    const int i(0);
    for (int j = 0; j <= p; j++)
    {
        Tmap[i][j] = 23.0;
    }
}

std::vector<std::vector<int>> puntexp(const std::vector<double>& data, const std::vector<double> differentials)
{
    std::vector< std::vector<int> > pd(2, std::vector<int> (4, 0));
    pd[0][0] = ceil(data[23]/differentials[0]) - 1; pd[1][0] = ceil(data[23]/differentials[0]);
    pd[0][1] = ceil(data[24]/differentials[1]) - 1; pd[1][1] = ceil(data[24]/differentials[1]);    
    pd[0][2] = ceil(data[25]/differentials[0]) - 1; pd[1][2] = ceil(data[25]/differentials[0]);
    pd[0][3] = ceil(data[26]/differentials[1]) - 1; pd[1][3] = ceil(data[26]/differentials[1]);
    return pd;
}