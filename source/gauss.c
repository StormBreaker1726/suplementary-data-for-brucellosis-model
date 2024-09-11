#include "../include/gauss.h"

real_cpu f(real_cpu x)
{
    return 0;
}

real_cpu gaussian_quadrature(real_cpu a, real_cpu b)
{
    real_cpu points[] = { 0,
                          (1.0f / 3.0f * sqrt(5.0f - 2.0f * sqrt(10.0f / 7.0f))),
                          -(1.0f / 3.0f * sqrt(5.0f - 2.0f * sqrt(10.0f / 7.0f))),
                          (1.0f / 3.0f * sqrt(5.0f + 2.0f * sqrt(10.0f / 7.0f))),
                          -(1.0f / 3.0f * sqrt(5.0f + 2.0f * sqrt(10.0f / 7.0f))) };

    real_cpu weights[] = { (128.0f / 225.0f),
                           ((322.0f + 13.0f * sqrt(70.0f)) / 900.0f),
                           ((322.0f + 13.0f * sqrt(70.0f)) / 900.0f),
                           ((322.0f - 13.0f * sqrt(70.0f)) / 900.0f),
                           ((322.0f - 13.0f * sqrt(70.0f)) / 900.0f) };

    real_cpu approx          = 0;
    real_cpu interval_change = (b - a) / 2.0f;
    real_cpu interval_center = (b + a) / 2.0f;

    for (size_t i = 0; i < 5; i++)
    {
        approx += weights[i] * f(interval_change * points[i] + interval_center);
    }

    approx *= interval_change;

    return approx;
}