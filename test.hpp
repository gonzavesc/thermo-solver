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

std::vector<double> readfiledat();
void setvalues(const std::vector<double>& C, std::vector<Material>& M);
void discret_points(const std::vector<double>& delta, const std::vector<double>& C, std::vector<std::vector<int>>& points);
void setcoeff(const std::vector<double>& diff, std::vector<Material>& M);