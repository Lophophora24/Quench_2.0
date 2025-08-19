#include "main.h"

void init()
{
    std::cout << "Starting init()....\n";

    fill();
    init_();                  // Thermal state with high temperature
    init_for_generating();

    //generate_initial_conditions_vacuum();   // Vacuum state

    //generate_initial_conditions_thermal();

    std::cout << "Exit init()....\n";
}

int main()
{
    clock_t time_start = clock();

    init();

    average();

    //calculate_phi_exact();

    calculate_observables();

    printf("time: %.2fs\n", (double)(clock() - time_start) / CLOCKS_PER_SEC);

    return 0;
}

