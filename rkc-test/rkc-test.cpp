#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#define real_cpu float

real_cpu C1     = 0.9f;
real_cpu C2     = 0.002f;
real_cpu C3     = 0.02f;
real_cpu C4     = 0.00086f;
real_cpu C5     = 0.00023f;
real_cpu C6     = 0.0082f;
real_cpu C7     = 0.37f;
real_cpu C8     = 1.0f;
real_cpu C9     = 3.2f;
real_cpu C10    = 0.36f;
real_cpu C11    = 0.04f;
real_cpu Bm_max = 50.0f;
real_cpu Bn0    = 50.0f;

real_cpu f_V(real_cpu V, real_cpu A)
{
    return C1 * V - C2 * V * A;
}

real_cpu f_Bn(real_cpu Bn, real_cpu V)
{
    return C3 * (Bn0 - Bn) - C4 * V * Bn;
}

real_cpu f_Ba(real_cpu Ba, real_cpu Bm, real_cpu Bn, real_cpu V)
{
    return C4 * Bn * V - C5 * Ba + C6 * Bm * V - C7 * Ba;
}

real_cpu f_Bm(real_cpu Bm, real_cpu Ba, real_cpu V)
{
    return C5 * Ba + C8 * Bm * (1 - (Bm / Bm_max)) - C6 * Bm * V;
}

real_cpu f_A(real_cpu A, real_cpu Ba, real_cpu Bm)
{
    return C9 * Ba + C10 * Bm - C11 * A;
}

void rkc(real_cpu *V, real_cpu *Bn, real_cpu *Ba, real_cpu *Bm, real_cpu *A, real_cpu dt, real_cpu t_max)
{
    unsigned n_points = (unsigned)(t_max / dt);

    real_cpu f1, f2, f3, f4;

    for (size_t i = 0; i < n_points - 1; i++)
    {
        // resolvendo vírus
        f1       = f_V(V[i], A[i]);
        f2       = f_V(V[i] + (dt / 2.0) * f1, A[i]);
        f3       = f_V(V[i] + (dt / 2.0) * f2, A[i]);
        f4       = f_V(V[i] + dt * f3, A[i]);
        V[i + 1] = V[i] + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4);

        // resolvendo célula B naive
        f1        = f_Bn(Bn[i], V[i]);
        f2        = f_Bn(Bn[i] + (dt / 2.0) * f1, V[i]);
        f3        = f_Bn(Bn[i] + (dt / 2.0) * f2, V[i]);
        f4        = f_Bn(Bn[i] + dt * f3, V[i]);
        Bn[i + 1] = Bn[i] + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4);

        // resolvendo célula B de memória
        f1        = f_Bm(Bm[i], Ba[i], V[i]);
        f2        = f_Bm(Bm[i] + (dt / 2.0) * f1, Ba[i], V[i]);
        f3        = f_Bm(Bm[i] + (dt / 2.0) * f2, Ba[i], V[i]);
        f4        = f_Bm(Bm[i] + dt * f3, Ba[i], V[i]);
        Bm[i + 1] = Bm[i] + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4);

        // resolvendo célula B de ativa
        f1        = f_Ba(Ba[i], Bm[i], Bn[i], V[i]);
        f2        = f_Ba(Ba[i] + (dt / 2.0) * f1, Bm[i], Bn[i], V[i]);
        f3        = f_Ba(Ba[i] + (dt / 2.0) * f2, Bm[i], Bn[i], V[i]);
        f4        = f_Ba(Ba[i] + dt * f3, Bm[i], Bn[i], V[i]);
        Ba[i + 1] = Ba[i] + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4);

        // resolvendo Anticorpo
        f1       = f_A(A[i], Ba[i], Bm[i]);
        f2       = f_A(A[i] + (dt / 2.0) * f1, Ba[i], Bm[i]);
        f3       = f_A(A[i] + (dt / 2.0) * f2, Ba[i], Bm[i]);
        f4       = f_A(A[i] + dt * f3, Ba[i], Bm[i]);
        A[i + 1] = A[i] + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4);
    }
}

int main(int argc, char **argv)
{
    real_cpu dt       = 0.0001;
    real_cpu t_max    = 120.0;
    unsigned n_points = (unsigned)(t_max / dt);

    std::cout << "Number of points in time = " << n_points << std::endl;

    real_cpu *V  = (real_cpu *)malloc(sizeof(real_cpu) * n_points);
    real_cpu *Bn = (real_cpu *)malloc(sizeof(real_cpu) * n_points);
    real_cpu *Ba = (real_cpu *)malloc(sizeof(real_cpu) * n_points);
    real_cpu *Bm = (real_cpu *)malloc(sizeof(real_cpu) * n_points);
    real_cpu *A  = (real_cpu *)malloc(sizeof(real_cpu) * n_points);

    for (size_t i = 0; i < n_points; i++)
    {
        V[i]  = 0.0f;
        Bn[i] = 0.0f;
        Ba[i] = 0.0f;
        Bm[i] = 0.0f;
        A[i]  = 0.0f;
    }

    V[0]  = 10.0f;
    Bn[0] = 350.0f;
    Ba[0] = 0.0f;
    Bm[0] = 0.0f;
    A[0]  = 0.0f;

    rkc(V, Bn, Ba, Bm, A, dt, t_max);

    std::ofstream V_file;
    std::ofstream Bn_file;
    std::ofstream Ba_file;
    std::ofstream Bm_file;
    std::ofstream A_file;

    V_file.open("cpp_results/V.txt");
    Bn_file.open("cpp_results/Bn.txt");
    Ba_file.open("cpp_results/Ba.txt");
    Bm_file.open("cpp_results/Bm.txt");
    A_file.open("cpp_results/A.txt");

    for (size_t i = 0; i < n_points; i++)
    {
        V_file << V[i] << std::endl;
        Bn_file << Bn[i] << std::endl;
        Ba_file << Ba[i] << std::endl;
        Bm_file << Bm[i] << std::endl;
        A_file << A[i] << std::endl;
    }

    free(V);
    free(Bn);
    free(Ba);
    free(Bm);
    free(A);

    return 0;
}