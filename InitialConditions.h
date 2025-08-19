#pragma once


#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>   
#include <chrono>
#include <fftw3.h>
#include <vector>
#include <complex>
#include <random>
#include <iostream>

#include "PARAMETERS.h"

extern double x[SIZE_X];

double eta(double);
