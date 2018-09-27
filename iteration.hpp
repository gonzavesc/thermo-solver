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

void iterate_gauss(const int& ini, const int& fini, const int& inj, const int& finj, std::vector<std::vector<double>>& Tmap,const std::vector<std::vector<double>>& Tmap_p , const Material& mat, double& rms);
void variable_temperature(std::vector<std::vector<double>>& Tmap, const double& Time, const int& ini, const int& fini, const int& inj,const int& finj);
void qflow(std::vector<std::vector<double>>& Tmap, const int& inj, const int& finj, const int& i, const std::vector<double>& diff,const Material mat);
void fluid_contact(std::vector<std::vector<double>>& Tmap, const int& ini, const int& fini, const int& j, const std::vector<double>& diff,const Material mat);
void gauss_seidel(std::vector<std::vector<double>>& Tmap, const std::vector<std::vector<double>>& Tmap_p, const std::vector<Material>& M,const std::vector<Material>& Mcont, const std::vector<Material>& Mcontd, const std::vector<double>& diff, const double& Time, const std::vector<std::vector<int>>& material_points, const int& N);