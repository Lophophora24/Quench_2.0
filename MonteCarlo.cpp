#include "MonteCarlo.h"

double phi_[SIZE_T][SIZE_X];                // Obnullit nado!!!

double energy_density_[SIZE_T][SIZE_X];
double energy_density_kin[SIZE_T][SIZE_X];   // кинетическая часть
double energy_density_inter[SIZE_T][SIZE_X]; // взаимодействующая часть
double energy_density_grad[SIZE_T][SIZE_X];  // часть с градиентом

double energy_[SIZE_T];
double energy_kin[SIZE_T];
double energy_inter[SIZE_T];
double energy_grad[SIZE_T];

double phi_0_0phi_t_0[SIZE_T];
double phi_0_0phi_0_x[SIZE_X];

void average()
{
	FILE* phi_sol;
	fopen_s(&phi_sol, "phi_sol.txt", "w+");

	for (int j = 0; j < SIZE_T; ++j) {
		energy_[j] = 0;
		energy_kin[j] = 0;
		energy_inter[j] = 0;
		energy_grad[j] = 0;
		for (int i = 0; i < SIZE_X; ++i) {
			phi_[j][i] = 0;
			energy_density_[j][i] = 0;
			energy_density_kin[j][i] = 0;
			energy_density_inter[j][i] = 0;
			energy_density_grad[j][i] = 0;
			
			phi_0_0phi_t_0[j] = 0;
			phi_0_0phi_0_x[i] = 0;
		}
	}

	for (int s = 1; s <= M; ++s) {
		solve_with_cond(s-1);

		// Fill the file with phi[t][x]
		std::cout << "solved eq # " << s << '\n';

		// убрать этот цикл когда будем писать в файл phi_sol.txt
		for (int j = 0; j < SIZE_T - 1; ++j) {

			for (int i = 0; i < SIZE_X; ++i) {

				phi_[j][i] += phi[j][i] / M;

				energy_density_kin[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2)) / M;

				energy_density_inter[j][i] += (0.25 * g * pow(phi[j][i], 4)) / M;

				energy_[j] += energy_density_[j][i] * h;
				energy_kin[j] += energy_density_kin[j][i] * h;
				energy_inter[j] += energy_density_inter[j][i] * h;
				energy_grad[j] += energy_density_grad[j][i] * h;

				
				


				// Check signs!!!
				if (i == SIZE_X - 1) {
					energy_density_[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2) + 0.5 * ((phi[j + 1][1] - phi[j + 1][i]) / h) * ((phi[j][1] - phi[j][i]) / h) + 0.5 * (F(j + 1, i) + F(j, i))) / M;

					energy_density_grad[j][i] += (0.5 * ((phi[j + 1][1] - phi[j + 1][i]) / h) * ((phi[j][1] - phi[j][i]) / h)) / M;
				
					
				}
				else {
					energy_density_[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2) + 0.5 * ((phi[j + 1][i + 1] - phi[j + 1][i]) / h) * ((phi[j][i + 1] - phi[j][i]) / h) + 0.5 * (F(j + 1, i) + F(j, i))) / M;

					

					energy_density_grad[j][i] += (0.5 * ((phi[j + 1][i + 1] - phi[j + 1][i]) / h) * ((phi[j][i + 1] - phi[j][i]) / h)) / M;

					
				}
			}

			phi_0_0phi_t_0[j] += (phi[0][0] * phi[j][0]) / M;

		}

		for (int i = 0; i < SIZE_X; ++i) {
			phi_0_0phi_0_x[i] += (phi[0][0] * phi[0][i]) / M;
		}

		/* Раскомментировать это когда будем писать в файл phi_sol.txt
		for (int j = 0; j < SIZE_T; ++j) {
			for (int i = 0; i < SIZE_X; ++i) {
				fprintf(phi_sol, "%20.3lf", phi[j][i]);
			}
			fprintf(phi_sol, "\n");
		}

		fprintf(phi_sol, "\n");
		*/
	}

	fclose(phi_sol);
}

void calculate_observables()
{
	FILE* phi_sol;
	fopen_s(&phi_sol, "phi_sol.txt", "r");

	for (int i = 0; i < M; ++i) {
		// somehow read phi[j][i] Mth array
		/* Раскомментировать это когда будем писать в файл phi_sol.txt
		for (int j = 0; j < SIZE_T; ++j) {
			for (int i = 0; i < SIZE_X; ++i) {
				fscanf_s(phi_sol, "%lf", &phi[j][i]);
			}
		}
		*/
		/*
		std::cout << "Read " << i << "th array\n";
		for (int j = 0; j < SIZE_T; ++j) {
			for (int i = 0; i < SIZE_X; ++i) {
				printf("%20.3lf", phi[j][i]);
			}
			printf("\n");
		}
		printf("\n");
		*/
		/* Раскомментировать это когда будем писать в файл phi_sol.txt
		for (int j = 0; j < SIZE_T-1; ++j) {
			for (int i = 0; i < SIZE_X; ++i) {

				phi_[j][i] += phi[j][i] / M;


				// Check signs!!!
				if (i == SIZE_X - 1) {
					energy_density_[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2) + 0.5 * ((phi[j + 1][1] - phi[j + 1][i]) / h) * ((phi[j][1] - phi[j][i]) / h) + 0.5 * (F(j + 1, i) + F(j, i))) / M;
				}
				else {
					energy_density_[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2) + 0.5 * ((phi[j + 1][i + 1] - phi[j + 1][i]) / h) * ((phi[j][i + 1] - phi[j][i]) / h) + 0.5 * (F(j + 1, i) + F(j, i))) / M;
				}
			}
		}
		*/
		
	}

	fclose(phi_sol);

	FILE* phi_aver;
	fopen_s(&phi_aver, "phi_aver.txt", "w+");

	FILE* energy_dens_aver;
	fopen_s(&energy_dens_aver, "energy_dens_aver.txt", "w+");
	FILE* energy_dens_kin_aver;
	fopen_s(&energy_dens_kin_aver, "energy_dens_kin_aver.txt", "w+");
	FILE* energy_dens_inter_aver;
	fopen_s(&energy_dens_inter_aver, "energy_dens_inter_aver.txt", "w+");
	FILE* energy_dens_grad_aver;
	fopen_s(&energy_dens_grad_aver, "energy_dens_grad_aver.txt", "w+");

	FILE* energy_aver;
	fopen_s(&energy_aver, "energy_aver.txt", "w+");

	FILE* phi_0_0phi_t_0_aver;
	fopen_s(&phi_0_0phi_t_0_aver, "phi_0_0phi_t_0_aver.txt", "w+");

	FILE* phi_0_0phi_0_x_aver;
	fopen_s(&phi_0_0phi_0_x_aver, "phi_0_0phi_0_x_aver.txt", "w+");

	for (int j = 0; j < SIZE_T; ++j) {
		for (int i = 0; i < SIZE_X; ++i) {
			fprintf(phi_aver, "%20.3lf", phi_[j][i]);
			fprintf(energy_dens_aver, "%20.3lf", energy_density_[j][i]);
			fprintf(energy_dens_kin_aver, "%20.3lf", energy_density_kin[j][i]);
			fprintf(energy_dens_inter_aver, "%20.3lf", energy_density_inter[j][i]);
			fprintf(energy_dens_grad_aver, "%20.3lf", energy_density_grad[j][i]);
		}
		fprintf(phi_aver, "\n");
		fprintf(energy_dens_aver, "\n");
		fprintf(energy_dens_kin_aver, "\n");
		fprintf(energy_dens_inter_aver, "\n");
		fprintf(energy_dens_grad_aver, "\n");

		fprintf(energy_aver, "%20.3lf %20.3lf %20.3lf %20.3lf\n", energy_[j], energy_kin[j], energy_inter[j], energy_grad[j]);

		fprintf(phi_0_0phi_t_0_aver, "%20.3lf\n", phi_0_0phi_t_0[j]);
	}

	for (int i = 0; i < SIZE_X; ++i) {
		fprintf(phi_0_0phi_0_x_aver, "%20.3lf\n", phi_0_0phi_0_x[i]);
	}

	fprintf(phi_aver, "\n");
	fprintf(energy_dens_aver, "\n");
	fprintf(energy_dens_kin_aver, "\n");
	fprintf(energy_dens_inter_aver, "\n");
	fprintf(energy_dens_grad_aver, "\n");
	fprintf(energy_aver, "\n");

	fclose(phi_aver);
	fclose(energy_dens_aver);
	fclose(energy_dens_kin_aver);
	fclose(energy_dens_inter_aver);
	fclose(energy_dens_grad_aver);
	fclose(energy_aver);
	fclose(phi_0_0phi_t_0_aver);
	fclose(phi_0_0phi_0_x_aver);

}