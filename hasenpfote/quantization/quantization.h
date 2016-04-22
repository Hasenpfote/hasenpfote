﻿/*!
* @file quantization.h
* @brief quantization helper.
* @author Hasenpfote
* @date 2016/04/23
*/
#pragma once
#include <cstdint>
//#include "../math/utility.h"

namespace hasenpfote{ namespace quantization{

// N-bit unsigned normalized value([0,1]) encoder.
template <unsigned N> std::uint16_t encode_unorm(float x)
{
    static_assert(((N) > 0) && ((N) <= 16), "out of range.");
    return static_cast<std::uint16_t>(x * ((1 << (N)) - 1) + 0.5f);
}

// N-bit unsigned normalized value([0,1]) decoder.
template <unsigned N> float decode_unorm(std::uint16_t x)
{
    static_assert(((N) > 0) && ((N) <= 16), "out of range.");
    return x / static_cast<float>((1 << (N)) - 1);
}

// N-bit signed normalized value([-1,+1]) encoder.
template <unsigned N> std::uint16_t encode_snorm(float x)
{
    static_assert(((N) > 1) && ((N) <= 16), "out of range.");
    return static_cast<std::uint16_t>((x < 0) | (encode_unorm<(N) - 1>(std::fabsf(x)) << 1));
}

// N-bit signed normalized value([-1,+1]) decoder.
template <unsigned N> float decode_snorm(std::uint16_t x)
{
    static_assert(((N) > 1) && ((N) <= 16), "out of range.");
    return decode_unorm<(N)- 1>(x >> 1) * ((x & 0x1) ? -1.0f : 1.0f);
}

}}