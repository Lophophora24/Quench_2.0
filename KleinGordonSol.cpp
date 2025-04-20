#include "KleinGordonSol.h"

double x[SIZE_X];
double time_[SIZE_T];

double phi[SIZE_T][SIZE_X];

void fill()
{
	for (int i = 0; i < SIZE_X; ++i)
		x[i] = x1 + i * h;
	for (int i = 0; i < SIZE_T; ++i)
		time_[i] = t1 + i * t;
}

void fill_phi()
{
	for (int i = 0; i < SIZE_X; ++i) {
		for (int j = 0; j < SIZE_T; ++j) {
			phi[j][i] = 0;
		}
		phi[0][i] = phi_0[i];
		phi[1][i] = phi[0][i] + pi_0[i] * t;
	}
}

double F(int j, int i)
{
	return m * m * pow(phi[j][i], 2) / 2 + g * pow(phi[j][i], 4) / 4;

    //return pow(m, 4) / (6 * g) * (cos(sqrt(6 * g) / m * phi[j][i]) - 1);
}

void Strauss_Vazquez()
{
	for (int j = 1; j < SIZE_T - 1; ++j) {
		for (int i = 0; i < SIZE_X; ++i) {

            phi[j + 1][i] = phi[j][i];
            double recent = phi[j + 1][i];

            if (i == 0) {
                do {
                    recent = phi[j + 1][i];
                    if (phi[j + 1][i] == phi[j - 1][i]) {
                        phi[j + 1][i] = r * r * phi[j][SIZE_X - 2] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                    }
                    else {
                        phi[j + 1][i] = r * r * phi[j][SIZE_X - 2] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                    }

                } while (fabs(recent - phi[j + 1][i]) > 1e-5);
            }
            else
            if (i == (SIZE_X - 1)) {
                do {
                    recent = phi[j + 1][i];
                    if (phi[j + 1][i] == phi[j - 1][i]) {
                       phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][1] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                    }
                    else {
                       phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][1] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                    }
                } while (fabs(recent - phi[j + 1][i]) > 1e-5);
            }
            else {
                do {
                    recent = phi[j + 1][i];
                    if (phi[j + 1][i] == phi[j - 1][i]) {
                        phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                    }
                    else {
                        phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                    }
                } while (fabs(recent - phi[j + 1][i]) > 1e-5);
            }
		}
	}
}

double eta(double x) {
	return 1. / (sqrt(2. * PI) * sigma) * exp(-(x - x_q) * (x - x_q) / (2. * sigma * sigma));
}

void solve_with_cond(int iter_num)
{
	generate_initial_conditions(iter_num);
	fill_phi();
	Strauss_Vazquez();
}

