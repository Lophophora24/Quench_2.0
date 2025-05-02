#include "main.h"

void init()
{
    std::cout << "Starting init()....\n";

    fill();
    init_();
    init_for_generating();

    std::cout << "Exit init()....\n";
}

int main()
{
    init();

    average();

    calculate_observables();

    return 0;
}

