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

void fill_phi(double **phi_local)
{
	for (int i = 0; i < SIZE_X; ++i) {
		for (int j = 0; j < SIZE_T; ++j) {
			phi_local[j][i] = 0;
		}
		phi_local[0][i] = phi_0[i];
		phi_local[1][i] = phi[0][i] + pi_0[i] * t;
	}
}

void print_phi(double** phi_local)
{
    for (int j = 0; j < SIZE_T; ++j) {
        for (int i = 0; i < SIZE_X; ++i) {
            std::cout << phi_local[j][i] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

double F(int j, int i, double** phi_local)
{
	return m * m * pow(phi_local[j][i], 2) / 2 + g * pow(phi_local[j][i], 4) / 4;

    //return pow(m, 4) / (6 * g) * (cos(sqrt(6 * g) / m * phi[j][i]) - 1);
}

void Strauss_Vazquez(double** phi_local)
{
    //print_phi(phi_local);
	for (int j = 1; j < SIZE_T - 1; ++j) {
		for (int i = 0; i < SIZE_X; ++i) {

            phi_local[j + 1][i] = phi_local[j][i];
            double recent = phi_local[j + 1][i];

            if (i == 0) {
                do {
                    recent = phi_local[j + 1][i];
                    if (phi_local[j + 1][i] == phi_local[j - 1][i]) {
                        phi_local[j + 1][i] = r * r * phi_local[j][SIZE_X - 2] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][i + 1] - phi_local[j - 1][i] - t * t * (m * m * phi_local[j + 1][i] + g * pow(phi_local[j + 1][i], 3));
                    }
                    else {
                        phi_local[j + 1][i] = r * r * phi_local[j][SIZE_X - 2] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][i + 1] - phi_local[j - 1][i] - t * t * (F(j + 1, i, phi_local) - F(j - 1, i, phi_local)) / (phi_local[j + 1][i] - phi_local[j - 1][i]);
                    }
                    //std::cout << phi_local[j + 1][i] << '\n';
                } while (fabs(recent - phi_local[j + 1][i]) > 1e-5);
            }
            else
            if (i == (SIZE_X - 1)) {
                do {
                    recent = phi_local[j + 1][i];
                    if (phi_local[j + 1][i] == phi_local[j - 1][i]) {
                       phi_local[j + 1][i] = r * r * phi_local[j][i - 1] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][1] - phi_local[j - 1][i] - t * t * (m * m * phi_local[j + 1][i] + g * pow(phi_local[j + 1][i], 3));
                    }
                    else {
                       phi_local[j + 1][i] = r * r * phi_local[j][i - 1] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][1] - phi_local[j - 1][i] - t * t * (F(j + 1, i, phi_local) - F(j - 1, i, phi_local)) / (phi_local[j + 1][i] - phi_local[j - 1][i]);
                    }
                    //std::cout << phi_local[j + 1][i] << '\n';
                } while (fabs(recent - phi_local[j + 1][i]) > 1e-5);
            }
            else {
                do {
                    recent = phi_local[j + 1][i];
                    if (phi_local[j + 1][i] == phi_local[j - 1][i]) {
                        phi_local[j + 1][i] = r * r * phi_local[j][i - 1] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][i + 1] - phi_local[j - 1][i] - t * t * (m * m * phi_local[j + 1][i] + g * pow(phi_local[j + 1][i], 3));
                    }
                    else {
                        phi_local[j + 1][i] = r * r * phi_local[j][i - 1] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][i + 1] - phi_local[j - 1][i] - t * t * (F(j + 1, i, phi_local) - F(j - 1, i, phi_local)) / (phi_local[j + 1][i] - phi_local[j - 1][i]);
                    }
                    //std::cout << phi_local[j + 1][i] << '\n';
                } while (fabs(recent - phi_local[j + 1][i]) > 1e-5);
            }
		}
	}
}

double eta(double x) {
	return 1. / (sqrt(2. * PI) * sigma) * exp(-(x - x_q) * (x - x_q) / (2. * sigma * sigma));
}

void solve_with_cond(int iter_num, double **phi_local)
{
    //std::cout << "starting solve_with_cond()....\n";
	
    generate_initial_conditions(iter_num);

    //std::cout << "starting fill_phi()\n";

	fill_phi(phi_local);

    //std::cout << "starting S-V()\n";
	
    Strauss_Vazquez(phi_local);
}



