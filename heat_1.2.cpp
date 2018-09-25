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
