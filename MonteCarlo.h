#pragma once

#include <iostream>
#include "PARAMETERS.h"

extern double x[SIZE_X];
extern double time_[SIZE_T];
extern double phi[SIZE_T][SIZE_X];
extern double phi_exact[SIZE_T][SIZE_X];

void solve_with_cond(int);

void print_phi(double**);

double F(int, int);