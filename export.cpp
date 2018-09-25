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

#ifndef INCLUDE_IT
    #include "iteration.hpp"
    #define INCLUDE_IT
#endif

#ifndef INCLUDE_EXP
    #include "export.hpp"
    #define INCLUDE_EXP
#endif

void exportarMatriu(const std::vector< std::vector<double> >& Tmap, const std::vector<std::vector<int>>& material_points)
{
    int i,jj;
    std::ofstream output;
    output.open("Matrix.dat");

    for (i = 0; i <= material_points[1][2]; i++)
    {
        for (jj = 0; jj <= material_points[0][2] - 1; jj++)
        {
             output << Tmap[i][jj] << ", " ;
        }
        output << Tmap[i][material_points[0][2]] << std::endl;
    }
    output.close();
}

void exportarTemp(const std::vector< std::vector<double> >& Tmap,const std::vector<std::vector<int>>& pd, const std::vector<double> differentials, const double& t, std::vector<double>& data)
{
    double T1, T2;
    std::ofstream output_T;
    output_T.open("Temp.dat", std::ios::out | std::ios::app);
    T1 = 1 / (differentials[0] * differentials[0]) * (Tmap[pd[0][1]][pd[0][0]] * (pd[1][0] * differentials[0] - data[23]) * (pd[1][1] * differentials[0] - data[24]) + 
    Tmap[pd[0][1]][pd[1][0]] * (data[23] - pd[0][0] * differentials[0]) * (pd[1][1] * differentials[0] - data[24]) + 
    Tmap[pd[1][1]][pd[0][0]] * (pd[1][0] * differentials[0] - data[23]) * (data[24] - pd[0][1] * differentials[0]) +
    Tmap[pd[1][1]][pd[1][0]] * (data[23] - pd[0][0] * differentials[0]) * (data[24] - pd[0][1] * differentials[0]));
    T2 = 1 / (differentials[0] * differentials[0]) * (Tmap[pd[0][3]][pd[0][2]] * (pd[1][2] * differentials[0] - data[25]) * (pd[1][3] * differentials[0] - data[26]) + 
    Tmap[pd[0][3]][pd[1][2]] * (data[25] - pd[0][2] * differentials[0]) * (pd[1][3] * differentials[0] - data[26]) + 
    Tmap[pd[1][3]][pd[0][2]] * (pd[1][2] * differentials[0] - data[25]) * (data[26] - pd[0][3] * differentials[0]) +
    Tmap[pd[1][3]][pd[1][2]] * (data[25] - pd[0][2] * differentials[0]) * (data[26] - pd[0][3] * differentials[0]));
    output_T << t <<", " << T1 << ", " << T2 << std::endl;
    output_T.close();
}
