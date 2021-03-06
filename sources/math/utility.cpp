﻿#include <cstdint>
#include <cmath>
#include <iostream>
#include <cstring>
#include "utility.h"

namespace hasenpfote{ namespace math{

bool almost_equals(float a, float b, std::uint32_t max_ulps)
{
    // Infinity check
    // FLT_MAX 近辺での比較を行わないよう単純比較に切り替える
    if(std::isinf(a) || std::isinf(b)){
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
        return a == b;
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }
    // NaN check
    // max_ulps が NaN との比較を行わないようにする
    // NaN - FLT_MAX - 1 でも代替可
    if(std::isnan(a) || std::isnan(b))
        return false;
    // Sign check
    //auto ia = *reinterpret_cast<std::int32_t*>(&a);
    //auto ib = *reinterpret_cast<std::int32_t*>(&b);
    std::int32_t ia, ib;
    std::memcpy(&ia, &a, sizeof(ia));
    std::memcpy(&ib, &b, sizeof(ib));

    if((ia & 0x80000000) != (ib & 0x80000000)){
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
        return a == b;
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }
    // 辞書順に並べ替え比較
    if(ia < 0)
        ia = 0x80000000 - ia;
    if(ib < 0)
        ib = 0x80000000 - ib;

    //std::cout << std::abs(ia - ib) << std::endl;

    return (static_cast<uint32_t>(std::abs(ia - ib)) <= max_ulps);
}

bool almost_equals(double a, double b, std::uint64_t max_ulps)
{
    // Infinity check
    // DBL_MAX 近辺での比較を行わないよう単純比較に切り替える
    if(std::isinf(a) || std::isinf(b)){
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
        return a == b;
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }
    // NaN check
    // max_ulps が NaN との比較を行わないようにする
    // NaN - DBL_MAX - 1 でも代替可
    if(std::isnan(a) || std::isnan(b))
        return false;
    // Sign check
    //auto ia = *reinterpret_cast<std::int64_t*>(&a);
    //auto ib = *reinterpret_cast<std::int64_t*>(&b);
    std::int64_t ia, ib;
    std::memcpy(&ia, &a, sizeof(ia));
    std::memcpy(&ib, &b, sizeof(ib));

    if((ia & 0x8000000000000000L) != (ib & 0x8000000000000000L)){
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
        return a == b;
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }
    // 辞書順に並べ替え比較
    if(ia < 0)
        ia = 0x8000000000000000L - ia;
    if(ib < 0)
        ib = 0x8000000000000000L - ib;

    std::cout << std::abs(ia - ib) << std::endl;

    return (static_cast<uint64_t>(std::abs(ia - ib)) <= max_ulps);
}

}}
