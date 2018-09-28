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

#ifndef INCLUDE_IT
    #include "iteration.hpp"
    #define INCLUDE_IT
#endif

void iterate_gauss(const int& ini, const int& fini, const int& inj, const int& finj, std::vector<std::vector<double>>& Tmap,
            const std::vector<std::vector<double>>& Tmap_p , const Material& mat, double& rms)
{
    int i, j;
    double prev;
    for (i = ini; i <= fini; i++)
    {
        for (j = inj; j <= finj; j++)
        {
            prev = Tmap[i][j];
            
            Tmap[i][j] = (mat.ae * Tmap[i][j+1] + mat.aw * Tmap[i][j-1] + mat.an * Tmap[i+1][j] + mat.as * Tmap[i-1][j] + mat.a0 * Tmap_p[i][j])
            / (mat.ap); 
            
            rms = rms + (prev - Tmap[i][j]) * (prev - Tmap[i][j]);
        }
    }
}

void variable_temperature(std::vector<std::vector<double>>& Tmap, const double& Time, const int& ini, const int& fini, const int& inj,
                          const int& finj)
{
    int i, j;
    for (i = ini; i <= fini; i++)
    {
        for (j = inj; j <= finj; j++)
        {
            Tmap[i][j] = 8 + 0.005 * Time;
        }
    }
}

void qflow(std::vector<std::vector<double>>& Tmap, const int& inj, const int& finj, const int& i, const std::vector<double>& diff,
           const Material mat)
{
    for (int j = inj; j <= finj; j++)
    {
        Tmap[i][j] = Tmap[i-1][j] + 60 * diff[1] / mat.k;
    }
}
void fluid_contact(std::vector<std::vector<double>>& Tmap, const int& ini, const int& fini, const int& j, const std::vector<double>& diff,
           const Material mat)
{
    for (int i = ini; i <= fini; i++)
    {
        Tmap[i][j] = (mat.k * Tmap[i][j+1] / diff[0] + 9.0 * 33.0) / (9.0 + mat.k / diff[0]);
    }
}

void gauss_seidel(std::vector<std::vector<double>>& Tmap, const std::vector<std::vector<double>>& Tmap_p, const std::vector<Material>& M,
                    const std::vector<Material>& Mcont, const std::vector<Material>& Mcontd, const std::vector<double>& diff, const double& Time, 
                    const std::vector<std::vector<int>>& material_points, const int& N)
{
    const double err (1.e-8);
    double rms;
    rms = 10.0;
    
    while( rms > err )
    {
        rms = 0;

        //boundary conditions
        //material 1
        
        //left
        fluid_contact(Tmap, 1, material_points[1][0], 0, diff, M[0]);
        
        //material 2
        //right
        variable_temperature(Tmap, Time, 1, material_points[1][1], material_points[0][2], material_points[0][2]);
        
        //material 3
        //top
        qflow(Tmap, 1, material_points[0][1], material_points[1][2], diff, M[2]);
        //left
        fluid_contact(Tmap, material_points[1][0] + 1, material_points[1][2] - 1, 0, diff, M[2]);
        
        //material 4
        //right
        variable_temperature(Tmap, Time, material_points[1][1] + 1, material_points[1][2], material_points[0][2], material_points[0][2]);
        
        //top
        
        qflow(Tmap, material_points[0][1] + 1, material_points[0][2] - 1, material_points[1][2], diff, M[3]);
        
        //material 1
        iterate_gauss(1, material_points[1][0] - 1, 1, material_points[0][0] - 1, Tmap, Tmap_p , M[0], rms);

        //material 2
        
        iterate_gauss(1, material_points[1][1] - 1, material_points[0][0] + 2, material_points[0][2] - 1, Tmap, Tmap_p , M[1], rms);
        
        //material 3
        
        
        iterate_gauss(material_points[1][0] + 2, material_points[1][2] - 1, 1, material_points[0][0] - 1, Tmap, Tmap_p , M[2], rms);
        //material 4
        iterate_gauss(material_points[1][1] + 2, material_points[1][2] - 1, material_points[0][0] + 2, material_points[0][2] - 1, Tmap, Tmap_p , M[3], rms);
        
        
        //boundary 1-2
        iterate_gauss(1, material_points[1][0] - 1, material_points[0][0], material_points[0][0], Tmap, Tmap_p , Mcont[0], rms);
        //boundary 1-2-3
        iterate_gauss(material_points[1][0], material_points[1][0], material_points[0][0], material_points[0][0], Tmap, Tmap_p , Mcontd[0], rms);
        
        
        
        //boundary 2-1
        iterate_gauss(1, material_points[1][0], material_points[0][0] + 1, material_points[0][0] + 1, Tmap, Tmap_p , Mcont[1], rms);
        
        
        
        //boundary 1-3
        iterate_gauss(material_points[1][0], material_points[1][0], 1, material_points[0][0] - 1, Tmap, Tmap_p , Mcont[2], rms);
        //boundary 3-1
        iterate_gauss(material_points[1][0] + 1, material_points[1][0] + 1, 1, material_points[0][0] - 1, Tmap, Tmap_p , Mcont[3], rms);
        //boundary 3-1-2
        iterate_gauss(material_points[1][0] + 1, material_points[1][0] + 1, material_points[0][0], material_points[0][0], Tmap, Tmap_p , Mcontd[2], rms);
        //boundary 2-3
        iterate_gauss(material_points[1][0] + 1, material_points[1][1] - 1, material_points[0][0] + 1, material_points[0][0] + 1, Tmap, Tmap_p , Mcont[4], rms);
    
        //boundary 3-2
        iterate_gauss(material_points[1][0] + 1, material_points[1][1], material_points[0][0], material_points[0][0], Tmap, Tmap_p , Mcont[5], rms);

        //boundary 2-3-4
        iterate_gauss(material_points[1][1], material_points[1][1], material_points[0][0] + 1, material_points[0][0] + 1, Tmap, Tmap_p , Mcontd[1], rms);
        //boundary 2-4
        iterate_gauss(material_points[1][1], material_points[1][1], material_points[0][0] + 2, material_points[0][2] - 1, Tmap, Tmap_p , Mcont[6], rms);
        //boundary 4-2
        iterate_gauss(material_points[1][1] + 1, material_points[1][1] + 1, material_points[0][0] + 2, material_points[0][2] - 1, Tmap, Tmap_p , Mcont[7], rms);
        //boundary 4-2-3
        iterate_gauss(material_points[1][1] + 1, material_points[1][1] + 1, material_points[0][0] + 1, material_points[0][0] + 1, Tmap, Tmap_p , Mcontd[3], rms);
        //boundary 3-4
        iterate_gauss(material_points[1][1] + 1, material_points[1][2] - 1, material_points[0][0], material_points[0][0], Tmap, Tmap_p , Mcont[8], rms);
        //boundary 4-3
        iterate_gauss(material_points[1][1] + 2, material_points[1][2] - 1, material_points[0][0] + 1, material_points[0][0] + 1, Tmap, Tmap_p , Mcont[8], rms);
        
        
    

        rms = sqrt(rms / N);    
          
    }

    


}