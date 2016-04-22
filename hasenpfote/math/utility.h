/**
* @file utility.h
* @brief math utility.
* @author Hasenpfote
* @date 2016/04/16
*/
#pragma once
#include <cstdint>

namespace hasenpfote{ namespace math{

/*!
 * a と b がほぼ等しいか
 * @param a
 * @param b
 * @param max_ulps unit in the last place.
 * @return a ≈ b
 */
bool almost_equals(float a, float b, std::uint32_t max_ulps);
bool almost_equals(double a, double b, std::uint64_t max_ulps);

/*!
 * sinc 関数のテイラー級数による近似.
 * <p>
 * sinc(x)<br>
 * \f$\frac{\sin{x}}{x} = 1 - \frac{1}{6}x^{2} + \frac{1}{120}x^{4} - \frac{1}{5040}x^{6}\f$
 * </p>
 * @param[in] x
 * @return
 */
float sinc(float x);

/*!
 * sinc 関数の逆数のテイラー級数による近似.
 * <p>
 * reciprocal of sinc(x)<br>
 * \f$\frac{x}{\sin{x}} = 1 + \frac{1}{6}x^{2} + \frac{7}{360}x^{4} + \frac{31}{15120}x^{6} + ･･･ \f$
 * </p>
 * @param x
 * @return
 */
float rcp_sinc(float x);

/*!
 * degrees to radians.
 * @param[in] angle an angle in degrees.
 * @return an angle in radians.
 */
float to_radians(float angle);

/*!
 * radians to degrees.
 * @param[in] angle an angle in radians.
 * @return an angle in degrees.
 */
float to_degrees(float angle);

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