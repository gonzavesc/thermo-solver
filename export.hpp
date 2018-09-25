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

void exportarMatriu(const std::vector< std::vector<double> >& Tmap, const std::vector<std::vector<int>>& material_points);
void exportarTemp(const std::vector< std::vector<double> >& Tmap,const std::vector<std::vector<int>>& pd, const std::vector<double> differentials, const double& t, std::vector<double>& data);