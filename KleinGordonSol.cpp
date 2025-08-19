#include "KleinGordonSol.h"

double x[SIZE_X];
double time_[SIZE_T];

double phi[SIZE_T][SIZE_X];

double phi_exact[SIZE_T][SIZE_X];

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

            phi_exact[j][i] = 0;
		}
		phi[0][i] = phi_0[i];
		phi[1][i] = phi[0][i] + pi_0[i] * t;
	}
}

void print_phi(double** phi_local)
{
    for (int j = 0; j < SIZE_T; ++j) {
        for (int i = 0; i < SIZE_X; ++i) {
            //std::cout << phi_local[j][i] << '     ';
            printf("%10.4lf", phi_local[j][i]);
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

double F(int j, int i)
{
	return m * m * pow(phi[j][i], 2) / 2 + g * pow(phi[j][i], 4) / 4;

    //return pow(m, 4) / (6 * g) * (cos(sqrt(6 * g) / m * phi[j][i]) - 1);
}
/*
using namespace Eigen;

void compute_F_and_J(VectorXd& u, VectorXd& u_prev, VectorXd& F, MatrixXd& J)
{
    for (int j = 0; j < SIZE_X; ++j) {
        int j_prev = (j - 1 + SIZE_X) % SIZE_X;
        int j_next = (j + 1) % SIZE_X;

        F(j) = (u[j] - 2 * u_prev[j] + u_prev[j]) / (t * t) -
            (u[j_next] - 2 * u[j] + u[j_prev]) / (h * h) +
            m * m * u[j] +
            g * u[j] * u[j] * u[j];

        J(j, j_prev) = -1. / (h * h);
        J(j, j) = 1. / (t * t) + 2. / (h * h) + m * m + 3 * g * u[j] * u[j];
        J(j, j_next) = -1. / (h * h);
    }
}

VectorXd Newton(VectorXd& u_curr, VectorXd& u_prev)
{
    VectorXd u_guess = u_curr;
    VectorXd F(SIZE_X);
    MatrixXd J(SIZE_X, SIZE_X);

    for (int iter = 0; iter < 20; ++iter) {
        J.setZero();
        compute_F_and_J(u_curr, u_prev, F, J);
        VectorXd delta = J.lu().solve(-F);
        u_guess += delta;

        if (delta.norm() < 1e-8) break;
    }
    return u_guess;
}

void Newton_solve(double** phi_local)
{
    VectorXd u_prev = VectorXd::Zero(SIZE_X);
    VectorXd u_curr = VectorXd::Zero(SIZE_X);

    for (int j = 0; j < SIZE_X; ++j) {
        u_prev[j] = phi_local[0][j];
        u_curr[j] = phi_local[1][j];
    }

    for (int i = 2; i < SIZE_T; ++i) {
        VectorXd u_next = Newton(u_curr, u_prev);
        u_prev = u_curr;
        u_curr = u_next;
        for (int j = 0; j < SIZE_X; ++j) {
            phi_local[i][j] = u_curr[j];
        }
    }

}*/

void Strauss_Vazquez()
{
    //print_phi(phi_local);
	for (int j = 1; j < SIZE_T - 1; ++j) {
		for (int i = 0; i < SIZE_X; ++i) {

            
            
            phi[j + 1][i] = phi[j][i];
            double recent = phi[j + 1][i];

            int i_prev = i - 1;
            int i_next = i + 1;

            if (i == 0) {
                i_prev = SIZE_X - 2;
            } 
            else if (i == SIZE_X - 1) {
                i_next = 1;
            }
            

            //double iter_array[20];

            for (int iter = 0; iter < 20; ++iter) {
                /*                
                if (std::isnan(phi_local[j+1][i])) {
                    for (int k = 0; k < iter; ++k) {
                        std::cout << "iter " << k << ": " << iter_array[k] << '\n';
                    }

                    getchar();
                }*/


                recent = phi[j + 1][i];
                
                if (phi[j + 1][i] == phi[j - 1][i]) {
                    phi[j + 1][i] = r * r * phi[j][i_prev] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i_next] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                }
                else {
                    phi[j + 1][i] = r * r * phi[j][i_prev] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i_next] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                }
                

                //phi_local[j + 1][i] = r * r * phi_local[j][i_prev] + 2 * (1 - r * r) * phi_local[j][i] + r * r * phi_local[j][i_next] - phi_local[j - 1][i] - t * t * 0.5 * (m * m * phi_local[j + 1][i] + g * pow(phi_local[j + 1][i], 3) + m * m * phi_local[j - 1][i] + g * pow(phi_local[j - 1][i], 3));

                //iter_array[iter] = phi_local[j + 1][i];



                if (fabs(recent - phi[j + 1][i]) < 1e-4) break;
            }
            
		}
	}
}

double eta(double x) {
	return 1. / (sqrt(2. * PI) * sigma) * exp(-(x - x_q) * (x - x_q) / (2. * sigma * sigma));
}

void solve_with_cond(int iter_num)
{
    //std::cout << "starting solve_with_cond()....\n";
	
    generate_initial_conditions(iter_num);

    //std::cout << "starting fill_phi()\n";

	fill_phi();

    //std::cout << "starting S-V()\n";
	
    Strauss_Vazquez();

    //Newton_solve(phi_local);
}

double theta_func(double x)
{
    if (x >= 0) return 1;
    else return 0;
}

double ret_green_func(double t_1, double t_2, double x_1, double x_2)
{
    double sum = 0;
    for (int i = -1000; i <= 1000; ++i) {
        if ((t_1 - t_2 - abs(x_1 - x_2 - i * L)) >= 0) {
            sum += _j0(m * sqrt((t_1 - t_2) * (t_1 - t_2) - (x_1 - x_2 - i * L) * (x_1 - x_2 - i * L)));
        }
        else {
            sum += 0;
        }
    }

    return -0.5 * sum;
}

void calculate_phi_exact()
{
    printf("Start calculate phi_exact....\n");
    for (int j = 0; j < SIZE_T; ++j) {
        if ((j % 100 == 0) || (j==(SIZE_T-1))) {
            for (int i = 0; i < SIZE_X; ++i) {
                for (int k = 0; k < SIZE_X; ++k) {
                    phi_exact[j][i] += al * h * eta(x[k]) * ret_green_func(time_[j], t1, x[i], x[k]);
                }
            }
        }

        printf("Finished iteration #%d\n", j);
    }
}



