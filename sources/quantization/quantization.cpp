#include "../assert.h"
#include "../fp.h"
#include "../math/constants.h"
#include "../math/utility.h"
#include "../math/quaternion.h"
#include "../math/vector3.h"
#include "quantization.h"

namespace hasenpfote{ namespace quantization{

/*!
 * 絶対値が最大の成分を取得.
 * @param[in] q
 * @return Quaternion::Component
 */
static math::Quaternion::Component get_component_by_abs_max(const math::Quaternion& q)
{
    auto component = math::Quaternion::Component::W;
    auto max = std::abs(q.GetW());
    auto abs = std::abs(q.GetX());
    if(abs > max){
        component = math::Quaternion::Component::X;
        max = abs;
    }
    abs = std::abs(q.GetY());
    if(abs > max){
        component = math::Quaternion::Component::Y;
        max = abs;
    }
    abs = std::abs(q.GetZ());
    if(abs > max){
        component = math::Quaternion::Component::Z;
    }
    return component;
}

std::uint32_t encode32_quat(const math::Quaternion& q)
{
    HASENPFOTE_ASSERT_MSG(math::almost_equals(1.0f, q.Norm(), 1), "Not an unit quaternion.");
    constexpr auto bit = 10;
    constexpr auto max = 1.0f / math::root_two<float>();    // + 1/sqrt(2)
    constexpr auto min = -max;                              // - 1/sqrt(2)
    std::uint32_t result = 0;
    const auto component = get_component_by_abs_max(q);
    const auto c = static_cast<std::underlying_type<math::Quaternion::Component>::type>(component);
    switch(component){
        case math::Quaternion::Component::W:
            {
                float x, y, z;
                if(q.GetW() < 0.0f){
                    x = math::remap(-q.GetX(), min, max, -1.0f, 1.0f);
                    y = math::remap(-q.GetY(), min, max, -1.0f, 1.0f);
                    z = math::remap(-q.GetZ(), min, max, -1.0f, 1.0f);
                }
                else{
                    x = math::remap(q.GetX(), min, max, -1.0f, 1.0f);
                    y = math::remap(q.GetY(), min, max, -1.0f, 1.0f);
                    z = math::remap(q.GetZ(), min, max, -1.0f, 1.0f);
                }
                result = (c << 30) | (encode_snorm<bit>(x) << 20) | (encode_snorm<bit>(y) << 10) | encode_snorm<bit>(z);
            }
            break;
        case math::Quaternion::Component::X:
            {
                float w, y, z;
                if(q.GetX() < 0.0f){
                    w = math::remap(-q.GetW(), min, max, -1.0f, 1.0f);
                    y = math::remap(-q.GetY(), min, max, -1.0f, 1.0f);
                    z = math::remap(-q.GetZ(), min, max, -1.0f, 1.0f);
                }
                else{
                    w = math::remap(q.GetW(), min, max, -1.0f, 1.0f);
                    y = math::remap(q.GetY(), min, max, -1.0f, 1.0f);
                    z = math::remap(q.GetZ(), min, max, -1.0f, 1.0f);
                }
                result = (c << 30) | (encode_snorm<bit>(w) << 20) | (encode_snorm<bit>(y) << 10) | encode_snorm<bit>(z);
            }
            break;
        case math::Quaternion::Component::Y:
            {
                float w, x, z;
                if(q.GetY() < 0.0f){
                    w = math::remap(-q.GetW(), min, max, -1.0f, 1.0f);
                    x = math::remap(-q.GetX(), min, max, -1.0f, 1.0f);
                    z = math::remap(-q.GetZ(), min, max, -1.0f, 1.0f);
                }
                else{
                    w = math::remap(q.GetW(), min, max, -1.0f, 1.0f);
                    x = math::remap(q.GetX(), min, max, -1.0f, 1.0f);
                    z = math::remap(q.GetZ(), min, max, -1.0f, 1.0f);
                }
                result = (c << 30) | (encode_snorm<bit>(w) << 20) | (encode_snorm<bit>(x) << 10) | encode_snorm<bit>(z);
            }
            break;
        case math::Quaternion::Component::Z:
            {
                float w, x, y;
                if(q.GetZ() < 0.0f){
                    w = math::remap(-q.GetW(), min, max, -1.0f, 1.0f);
                    x = math::remap(-q.GetX(), min, max, -1.0f, 1.0f);
                    y = math::remap(-q.GetY(), min, max, -1.0f, 1.0f);
                }
                else{
                    w = math::remap(q.GetW(), min, max, -1.0f, 1.0f);
                    x = math::remap(q.GetX(), min, max, -1.0f, 1.0f);
                    y = math::remap(q.GetY(), min, max, -1.0f, 1.0f);
                }
                result = (c << 30) | (encode_snorm<bit>(w) << 20) | (encode_snorm<bit>(x) << 10) | encode_snorm<bit>(y);
            }
            break;
        default:
            HASENPFOTE_ASSERT_MSG(false, "Unknown type.");
            break;
    }
    return result;
}

math::Quaternion decode32_quat(std::uint32_t q)
{
    constexpr auto bit = 10;
    constexpr auto max = 1.0f / math::root_two<float>();    // + 1/sqrt(2)
    constexpr auto min = -max;                              // - 1/sqrt(2)
    float w, x, y, z;
    const auto component = static_cast<math::Quaternion::Component>(q >> 30);
    switch(component){
        case math::Quaternion::Component::W:
            {
                x = decode_snorm<bit>((q >> 20) & 0x3FF);
                y = decode_snorm<bit>((q >> 10) & 0x3FF);
                z = decode_snorm<bit>(q & 0x3FF);
                x = math::remap(x, -1.0f, 1.0f, min, max);
                y = math::remap(y, -1.0f, 1.0f, min, max);
                z = math::remap(z, -1.0f, 1.0f, min, max);
                w = std::sqrt(1.0f - x * x - y * y - z * z);
            }
            break;
        case math::Quaternion::Component::X:
            {
                w = decode_snorm<bit>((q >> 20) & 0x3FF);
                y = decode_snorm<bit>((q >> 10) & 0x3FF);
                z = decode_snorm<bit>(q & 0x3FF);
                w = math::remap(w, -1.0f, 1.0f, min, max);
                y = math::remap(y, -1.0f, 1.0f, min, max);
                z = math::remap(z, -1.0f, 1.0f, min, max);
                x = std::sqrt(1.0f - w * w - y * y - z * z);
            }
            break;
        case math::Quaternion::Component::Y:
            {
                w = decode_snorm<bit>((q >> 20) & 0x3FF);
                x = decode_snorm<bit>((q >> 10) & 0x3FF);
                z = decode_snorm<bit>(q & 0x3FF);
                w = math::remap(w, -1.0f, 1.0f, min, max);
                x = math::remap(x, -1.0f, 1.0f, min, max);
                z = math::remap(z, -1.0f, 1.0f, min, max);
                y = std::sqrt(1.0f - w * w - x * x - z * z);
            }
            break;
        case math::Quaternion::Component::Z:
            {
                w = decode_snorm<bit>((q >> 20) & 0x3FF);
                x = decode_snorm<bit>((q >> 10) & 0x3FF);
                y = decode_snorm<bit>(q & 0x3FF);
                w = math::remap(w, -1.0f, 1.0f, min, max);
                x = math::remap(x, -1.0f, 1.0f, min, max);
                y = math::remap(y, -1.0f, 1.0f, min, max);
                z = std::sqrt(1.0f - w * w - x * x - y * y);
            }
            break;
        default:
            HASENPFOTE_ASSERT_MSG(false, "Unknown type.");
            break;
    }
    return math::Quaternion(w, x, y, z);
}

std::uint64_t encode161616_vec(const math::Vector3& v)
{
    return (static_cast<std::uint64_t>(float_to_half(v.GetX())) << 32)
        | (static_cast<std::uint64_t>(float_to_half(v.GetY())) << 16)
        | static_cast<std::uint64_t>(float_to_half(v.GetZ()));
}

math::Vector3 decode161616_vec(std::uint64_t v)
{
    return math::Vector3(
        half_to_float((v >> 32) & 0xFFFF),
        half_to_float((v >> 16) & 0xFFFF),
        half_to_float(v & 0xFFFF));
}

}}