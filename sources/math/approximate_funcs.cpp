#include "approximate_funcs.h"

namespace hasenpfote{ namespace math{

float sinc(float x)
{
    const float x2 = x * x;
    const float x4 = x2 * x2;
    return 1.0f - x2 / 6.0f + x4 / 120.0f - (x2 * x4) / 5040.0f;
}

float rcp_sinc(float x)
{
    const float x2 = x * x;
    const float x4 = x2 * x2;
    return 1.0f + x2 / 6.0f + 7.0f * x4 / 360.0f + 31.0f * (x2 * x4) / 15120.0f;
}

}}