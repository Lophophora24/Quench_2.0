#include "main.h"

void init()
{
    fill();
    init_();
    init_for_generating();
}

int main()
{
    init();

    average();

    calculate_observables();

    return 0;
}

