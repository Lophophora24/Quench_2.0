#include "MonteCarlo.h"
#include "omp.h"

double phi_[SIZE_T][SIZE_X];                 // <φ(t,x)>

											 // ε(t,x) = ε_kin(t,x) + ε_grad(t,x) + ε_inter(t,x)
double energy_density_[SIZE_T][SIZE_X];      // <ε(t,x)> 
double energy_density_kin[SIZE_T][SIZE_X];   // <ε_kin(t,x)>
double energy_density_inter[SIZE_T][SIZE_X]; // <ε_inter(t,x)>
double energy_density_grad[SIZE_T][SIZE_X];  // <ε_grad(t,x)>

											 // E = E_kin + E_grad + E_inter
double energy_[SIZE_T];						 // E
double energy_kin[SIZE_T];                   // E_kin
double energy_inter[SIZE_T];                 // E_inter
double energy_grad[SIZE_T];                  // E_grad
						
											 // Correlation functions
double phi_0_0phi_t_0[SIZE_T];               // <φ(0,0)φ(t,0)>
double phi_0_0phi_0_x[SIZE_X];               // <φ(0,0)φ(0,x)>



void average()
{
	std::cout << "Starting average()....\n";
	// Obnulenie vseh massivov //
	for (int j = 0; j < SIZE_T; ++j) {
		energy_[j]      = 0;
		energy_kin[j]   = 0;
		energy_inter[j] = 0;
		energy_grad[j]  = 0;
		for (int i = 0; i < SIZE_X; ++i) {
			phi_[j][i]                 = 0;
			energy_density_[j][i]      = 0;

			energy_density_kin[j][i]   = 0;
			energy_density_inter[j][i] = 0;
			energy_density_grad[j][i]  = 0;
			
			phi_0_0phi_t_0[j]          = 0;
			phi_0_0phi_0_x[i]          = 0;

			
		}
	}

	//std::cout << "Starting parallel section....\n";
	//#pragma omp parallel reduction (+:phi_)
	 
		//#pragma omp parallel for		
		for (int s = 1; s <= M; ++s) {
			
			//double phi[SIZE_T][SIZE_X];
			/*
			double** phi;

			phi = (double**)malloc(SIZE_T * sizeof(double*));
			//std::cout << phi << '\n';
			for (int j = 0; j < SIZE_T; ++j) {
				phi[j] = (double*)malloc(SIZE_X * sizeof(double));
				for (int i = 0; i < SIZE_X; ++i) {
					phi[j][i] = 0;
				}
			}
			*/
			solve_with_cond(s - 1);

			// Fill the file with phi[t][x]
			std::cout << "solved eq # " << s << '\n';
			/*
			double en[SIZE_T];
			for (int j = 0; j < SIZE_T - 1; ++j) {
				en[j] = 0;
			}
			for (int j = 0; j < SIZE_T - 1; ++j) {
				for (int i = 0; i < SIZE_X; ++i) {
					int i_prev = i - 1;
					int i_next = i + 1;

					if (i == 0) {
						i_prev = SIZE_X - 2;
					}
					else if (i == SIZE_X - 1) {
						i_next = 1;
					}
					en[j] += h * (0.5 * (phi[j + 1][i] - phi[j][i]) / t + 0.5 * (phi[j + 1][i_next] - phi[j + 1][i]) / h * (phi[j][i_next] - phi[j][i]) / h + 0.5 * (F(j + 1, i) + F(j, i)));
				}
			}

			for (int j = 0; j < SIZE_T; ++j) {
				printf("%10.4lf\n", en[j]);
			}*/

			//print_phi(phi);
			
			//std::cout << omp_get_thread_num() << '\n';
			//std::cout << omp_get_max_threads() << '\n';

			//#pragma omp parallel for
			for (int j = 0; j < SIZE_T - 1; ++j) {
				for (int i = 0; i < SIZE_X; ++i) {
					int i_prev = i - 1;
					int i_next = i + 1;

					if (i == 0) {
						i_prev = SIZE_X - 2;
					}
					else if (i == SIZE_X - 1) {
						i_next = 1;
					}

					//#pragma omp atomic 
					phi_[j][i] += phi[j][i] / M;

					//energy_density_kin[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2)) / M;

					//energy_density_inter[j][i] += (0.25 * g * pow(phi[j][i], 4)) / M;

					energy_density_[j][i] += (0.5 * (phi[j + 1][i] - phi[j][i]) / t * (phi[j + 1][i] - phi[j][i]) / t + 0.5 * (phi[j + 1][i_next] - phi[j + 1][i]) / h * (phi[j][i_next] - phi[j][i]) / h + 0.5 * (F(j + 1, i) + F(j, i))) / M;

					energy_density_kin[j][i] += (0.5 * (phi[j + 1][i] - phi[j][i]) / t * (phi[j + 1][i] - phi[j][i]) / t) / M;

					energy_density_grad[j][i] += (0.5 * (phi[j + 1][i_next] - phi[j + 1][i]) / h * (phi[j][i_next] - phi[j][i]) / h) / M;

					energy_density_inter[j][i] += (0.5 * (F(j + 1, i) + F(j, i))) / M;

					//energy_density_[j][i] += (energy_density_kin[j][i] + energy_density_grad[j][i] + energy_density_inter[j][i]);

					//energy_[j] += energy_density_[j][i] * h;
					//energy_kin[j] += energy_density_kin[j][i] * h;
					//energy_inter[j] += energy_density_inter[j][i] * h;
					//energy_grad[j] += energy_density_grad[j][i] * h;

					/*
					if (i == SIZE_X - 1) {
						energy_density_[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2) + 0.5 * ((phi[j + 1][1] - phi[j + 1][i]) / h) * ((phi[j][1] - phi[j][i]) / h) + 0.5 * (F(j + 1, i, phi) + F(j, i, phi))) / M;
						//energy_density_grad[j][i] += (0.5 * ((phi[j + 1][1] - phi[j + 1][i]) / h) * ((phi[j][1] - phi[j][i]) / h)) / M;
					}
					else {
						energy_density_[j][i] += (0.5 * pow(((phi[j + 1][i] - phi[j][i]) / t), 2) + 0.5 * ((phi[j + 1][i + 1] - phi[j + 1][i]) / h) * ((phi[j][i + 1] - phi[j][i]) / h) + 0.5 * (F(j + 1, i, phi) + F(j, i, phi))) / M;
						//energy_density_grad[j][i] += (0.5 * ((phi[j + 1][i + 1] - phi[j + 1][i]) / h) * ((phi[j][i + 1] - phi[j][i]) / h)) / M;
					}
					*/
					
				}
				phi_0_0phi_t_0[j] += (phi[0][0] * phi[j][0]) / M;
			}

			for (int j = 0; j < SIZE_T - 1; ++j) {
				for (int i = 0; i < SIZE_X - 1; ++i) {
					
					energy_[j] += energy_density_[j][i] * h;

					energy_kin[j] += energy_density_kin[j][i] * h;
					energy_grad[j] += energy_density_grad[j][i] * h;
					energy_inter[j] += energy_density_inter[j][i] * h;

					//energy_[j] += (energy_kin[j] + energy_grad[j] + energy_inter[j]);
				}
			}

			for (int i = 0; i < SIZE_X; ++i) {
				phi_0_0phi_0_x[i] += (phi[0][0] * phi[0][i]) / M;
			} 


			for (int i = 0; i < SIZE_X; ++i) {
				phi_[SIZE_T - 1][i] = phi_[SIZE_T - 2][i];
				energy_density_[SIZE_T - 1][i] = energy_density_[SIZE_T - 2][i];
				energy_density_kin[SIZE_T - 1][i] = energy_density_kin[SIZE_T - 2][i];
				energy_density_grad[SIZE_T - 1][i] = energy_density_grad[SIZE_T - 2][i];
				energy_density_inter[SIZE_T - 1][i] = energy_density_inter[SIZE_T - 2][i];
			}
			energy_[SIZE_T - 1] = energy_[SIZE_T - 2];
			energy_kin[SIZE_T - 1] = energy_kin[SIZE_T - 2];
			energy_grad[SIZE_T - 1] = energy_grad[SIZE_T - 2];
			energy_inter[SIZE_T - 1] = energy_inter[SIZE_T - 2];

			phi_0_0phi_t_0[SIZE_T - 1] = phi_0_0phi_t_0[SIZE_T - 2];
			

		}
	
}

void calculate_observables()
{
	FILE* phi_aver;
	fopen_s(&phi_aver, "phi_aver.txt", "w+");

	FILE* phi_exact_;
	fopen_s(&phi_exact_, "phi_exact.txt", "w+");

	FILE* energy_dens_aver;
	fopen_s(&energy_dens_aver, "energy_dens_aver.txt", "w+");
	//FILE* energy_dens_kin_aver;
	//fopen_s(&energy_dens_kin_aver, "energy_dens_kin_aver.txt", "w+");
	//FILE* energy_dens_inter_aver;
	//fopen_s(&energy_dens_inter_aver, "energy_dens_inter_aver.txt", "w+");
	//FILE* energy_dens_grad_aver;
	//fopen_s(&energy_dens_grad_aver, "energy_dens_grad_aver.txt", "w+");

	FILE* energy_aver;
	fopen_s(&energy_aver, "energy_aver.txt", "w+");

	FILE* phi_0_0phi_t_0_aver;
	fopen_s(&phi_0_0phi_t_0_aver, "phi_0_0phi_t_0_aver.txt", "w+");

	FILE* phi_0_0phi_0_x_aver;
	fopen_s(&phi_0_0phi_0_x_aver, "phi_0_0phi_0_x_aver.txt", "w+");

	fprintf(energy_dens_aver, "%20.3lf", " ");
	fprintf(phi_aver, "%20.3lf", " ");
	fprintf(phi_exact_, "%20.3lf", " ");
	for (int j = 0; j < SIZE_T; ++j) {
		
		fprintf(energy_dens_aver, "%20.3lf", time_[j]);
		fprintf(phi_aver, "%20.3lf", time_[j]);
		fprintf(phi_exact_, "%20.3lf", time_[j]);
	}
	fprintf(energy_dens_aver, "\n");
	fprintf(phi_aver, "\n");
	fprintf(phi_exact_, "\n");

	for (int i = 0; i < SIZE_X; ++i) {
		fprintf(phi_aver, "%10.3lf", x[i]);
		fprintf(phi_exact_, "%10.3lf", x[i]);
		fprintf(energy_dens_aver, "%10.3lf", x[i]);
		for (int j = 0; j < SIZE_T; ++j) {
			fprintf(phi_aver, "%20.3lf", phi_[j][i]);
			fprintf(phi_exact_, "%20.3lf", phi_exact[j][i]);
			fprintf(energy_dens_aver, "%20.3lf", energy_density_[j][i]);
			//fprintf(energy_dens_kin_aver, "%20.3lf", energy_density_kin[j][i]);
			//fprintf(energy_dens_inter_aver, "%20.3lf", energy_density_inter[j][i]);
			//fprintf(energy_dens_grad_aver, "%20.3lf", energy_density_grad[j][i]);
		}
		fprintf(phi_aver, "\n");
		fprintf(phi_exact_, "\n");
		fprintf(energy_dens_aver, "\n");

		//fprintf(energy_dens_kin_aver, "\n");
		//fprintf(energy_dens_inter_aver, "\n");
		//fprintf(energy_dens_grad_aver, "\n");

		//fprintf(energy_aver, "%20.3lf %20.3lf %20.3lf %20.3lf\n", energy_[j], energy_kin[j], energy_inter[j], energy_grad[j]);

		//fprintf(phi_0_0phi_t_0_aver, "%20.3lf\n", phi_0_0phi_t_0[j]);
	}



	for (int j = 0; j < SIZE_T; ++j) {
		fprintf(energy_aver, "%10.3lf%20.3lf%20.3lf%20.3lf%20.3lf\n", time_[j], energy_[j], energy_kin[j], energy_grad[j], energy_inter[j]);
		fprintf(phi_0_0phi_t_0_aver, "%10.3lf%20.3lf\n", time_[j], phi_0_0phi_t_0[j]);

	}

	for (int i = 0; i < SIZE_X; ++i) {
		fprintf(phi_0_0phi_0_x_aver, "%10.3lf%20.3lf\n", x[i], phi_0_0phi_0_x[i]);
	}

	//fprintf(phi_aver, "\n");
	//fprintf(energy_dens_aver, "\n");
	//fprintf(energy_dens_kin_aver, "\n");
	//fprintf(energy_dens_inter_aver, "\n");
	//fprintf(energy_dens_grad_aver, "\n");
	//fprintf(energy_aver, "\n");

	fclose(phi_aver);
	fclose(phi_exact_);
	fclose(energy_dens_aver);
	//fclose(energy_dens_kin_aver);
	//fclose(energy_dens_inter_aver);
	//fclose(energy_dens_grad_aver);
	fclose(energy_aver);
	fclose(phi_0_0phi_t_0_aver);
	fclose(phi_0_0phi_0_x_aver);
}