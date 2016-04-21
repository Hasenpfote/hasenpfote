/**
* @file mathutil.h
* @brief math utility class.
* @author Hasenpfote
* @date 2016/04/16
*/
#pragma once
#include <iostream>
#include <cstdint>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

namespace hasenpfote{ namespace math{

class MathUtil final
{
private:
    MathUtil() = delete;
    ~MathUtil() = delete;

public:

    /*!
     * a と b がほぼ等しいか
     * @param a
     * @param b
     * @param max_ulps unit in the last place.
     * @return a ≈ b
     */
    static bool AlmostEquals(float a, float b, uint32_t max_ulps);

    /*!
     * sinc 関数のテイラー級数による近似.
     * <p>
     * sinc(x)<br>
     * \f$\frac{\sin{x}}{x} = 1 - \frac{1}{6}x^{2} + \frac{1}{120}x^{4} - \frac{1}{5040}x^{6}\f$
     * </p>
     * @param[in] x
     * @return
     */
    static float Sinc(float x);

    /*!
     * sinc 関数の逆数のテイラー級数による近似.
     * <p>
     * reciprocal of sinc(x)<br>
     * \f$\frac{x}{\sin{x}} = 1 + \frac{1}{6}x^{2} + \frac{7}{360}x^{4} + \frac{31}{15120}x^{6} + ･･･ \f$
     * </p>
     * @param x
     * @return
     */
    static float ReciprocalOfSinc(float x);

    static float DegreesToRadians(float angle);
    static float RadiansToDegrees(float angle);
};

inline float MathUtil::Sinc(float x)
{
    const float x2 = x * x;
    const float x4 = x2 * x2;
    return 1.0f - x2 / 6.0f + x4 / 120.0f - (x2 * x4) / 5040.0f;
}

inline float MathUtil::ReciprocalOfSinc(float x)
{
    const float x2 = x * x;
    const float x4 = x2 * x2;
    return 1.0f + x2 / 6.0f + 7.0f * x4 / 360.0f + 31.0f * (x2 * x4) / 15120.0f;
}

inline float MathUtil::DegreesToRadians(float angle)
{
    return angle * static_cast<float>(M_PI) / 180.0f;
}

inline float MathUtil::RadiansToDegrees(float angle)
{
    return angle * 180.0f / static_cast<float>(M_PI);
}

template <typename T>
T clamp(T value, T min, T max)
{
    if(value < min)
        return min;
    if(value > max)
        return max;
    return value;
};

template <typename T>
bool contains_closed(T value, T lower, T upper)
{
    if(value <= lower || value >= upper)
        return false;
    return true;
};

template <typename T>
bool contains_open(T value, T lower, T upper)
{
    if(value < lower || value > upper)
        return false;
    return true;
};


}}