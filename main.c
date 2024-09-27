#include "include/defines.h"
#include "include/inate_solver.h"

int main(int argc, char **argv)
{
    unsigned omp_threads_qtd = 1;

    unsigned seed = -1041698816;

    if (argc > 2)
    {
        omp_threads_qtd = atoi(argv[1]);
    }
    else if (argc > 1)
    {
        omp_threads_qtd = atoi(argv[1]);
    }

    printf("Number of threads: %d\n", omp_threads_qtd);
    printf("Seed: %d\n", seed);

    FILE *results_info = fopen("results/info.txt", "w");

    fprintf(results_info, "seed used in this experiment = %d\n", seed);
    solve_model(60, 100, 100, omp_threads_qtd, seed);

    return 0;
}
