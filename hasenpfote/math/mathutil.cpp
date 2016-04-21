#include <limits>
#include "mathutil.h"

namespace hasenpfote{ namespace math{

bool MathUtil::AlmostEquals(float a, float b, uint32_t max_ulps)
{
    // Infinity check
    // FLT_MAX 近辺での比較を行わないよう単純比較に切り替える
    if(std::isinf(a) || std::isinf(b))
        return a == b;
    // NaN check
    // max_ulps が NaN との比較を行わないようにする
    // NaN - FLT_MAX - 1 でも代替可
    if(std::isnan(a) || std::isnan(b))
        return false;
    // Sign check
    int32_t ia = *reinterpret_cast<uint32_t*>(&a);
    int32_t ib = *reinterpret_cast<uint32_t*>(&b);
    if((ia & 0x80000000) != (ib & 0x80000000))
        return a == b;
    // 辞書順に並べ替え比較
    if(ia < 0)
        ia = 0x80000000 - ia;
    if(ib < 0)
        ib = 0x80000000 - ib;

    //std::cout << std::abs(ia - ib) << std::endl;

    return (static_cast<uint32_t>(std::abs(ia - ib)) <= max_ulps);
}

}}