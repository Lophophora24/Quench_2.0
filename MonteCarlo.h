#pragma once

#include <iostream>
#include "PARAMETERS.h"

extern double x[SIZE_X];
extern double phi[SIZE_T][SIZE_X];

void solve_with_cond(int, double**);

double F(int, int, double**);