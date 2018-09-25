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
    const double err (1.e-4);
    double rms;
    rms = 10.0;
    double a (0);
    
    while( rms > err )
    {
        a = a-1;
        rms = 0;
        
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



}

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

int main()
{
    double impl, runtime;
    double Time;
    int n(0), N, i, j;
    std::vector<double> data;
    std::vector<double> differentials;
    Material mat_1, mat_2, mat_3, mat_4;
    Material contact_12, contact_21, contact_13, contact_31, contact_23, contact_32, contact_24, contact_42, contact_34, contact_43;
    Material contact_123, contact_234, contact_312, contact_423;
    std::vector<Material> M;
    std::vector<Material> Mcont, Mcontd;
    std::vector<std::vector<int>> material_points(2, std::vector<int>(3,0));
    std::vector<std::vector<int>> punts_discr;

    data = readfiledat();
    impl = data[22]; runtime = data[21];
    differentials.push_back(data[18]); differentials.push_back(data[19]); differentials.push_back(data[20]);
    punts_discr = puntexp(data, differentials);

    M.push_back(mat_1); M.push_back(mat_2); M.push_back(mat_3); M.push_back(mat_4);

    setvalues(data, M);
    discret_points(differentials, data, material_points);
    
    setcoeff(differentials, M);

    setcontact(differentials, contact_12, M[0], M[1], 0, 0); Mcont.push_back(contact_12);
    setcontact(differentials, contact_21, M[1], M[0], 0, 1); Mcont.push_back(contact_21);
    setcontact(differentials, contact_13, M[0], M[2], 1, 0); Mcont.push_back(contact_13);
    setcontact(differentials, contact_31, M[2], M[0], 1, 1); Mcont.push_back(contact_31);
    setcontact(differentials, contact_23, M[1], M[2], 0, 1); Mcont.push_back(contact_23);
    setcontact(differentials, contact_32, M[2], M[1], 0, 0); Mcont.push_back(contact_32);
    setcontact(differentials, contact_24, M[1], M[3], 1, 0); Mcont.push_back(contact_24);
    setcontact(differentials, contact_42, M[3], M[1], 1, 1); Mcont.push_back(contact_42);
    setcontact(differentials, contact_34, M[2], M[3], 0, 0); Mcont.push_back(contact_34);
    setcontact(differentials, contact_43, M[3], M[2], 0, 1); Mcont.push_back(contact_43);

    setcontact(differentials, contact_123, contact_12, M[2], 1, 0); Mcontd.push_back(contact_123);
    setcontact(differentials, contact_234, contact_23, M[3], 1, 0); Mcontd.push_back(contact_234);
    setcontact(differentials, contact_312, contact_31, M[1], 0, 0); Mcontd.push_back(contact_312);
    setcontact(differentials, contact_423, contact_42, M[2], 0, 1); Mcontd.push_back(contact_423);

    std::vector<std::vector<double>> Tmap(material_points[1][2] + 1, std::vector<double>(material_points[0][2] + 1,8));
    std::vector<std::vector<double>> Tmap_p(material_points[1][2] + 1, std::vector<double>(material_points[0][2] + 1,8));
    initializeT(Tmap, material_points[0][2]);

    N = (material_points[1][2] - 1) * (material_points[0][2] - 1);


    Time = n * differentials[3];
    if (impl == 0.0)
    {

    }
    if (impl == 0.5)
    {

    }
    if (impl == 1)
    {
        while(Time < runtime)
        {
            n++;
            Time = differentials[2] * n;
            for (i = 0; i <= material_points[1][2]; i++)
            {
                for (j = 0; j <= material_points[0][2]; j++)
                {
                    Tmap_p[i][j] = Tmap[i][j];
                }
            }
            exportarTemp(Tmap, punts_discr, differentials, Time, data);
            gauss_seidel(Tmap, Tmap_p, M, Mcont, Mcontd, differentials, Time, material_points, N);
        }
        exportarMatriu(Tmap, material_points);
    }
    return 0;
}
