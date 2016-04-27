/*!
* @file fp.h
* @brief Floating point number.
* @author Hasenpfote
* @date 2016/04/27
*/
#pragma once
#include <cstdint>
#include <ostream>

namespace hasenpfote{

union FP32;

union FP16
{
    std::uint16_t u;
    struct
    {
        std::uint16_t mantissa : 10;
        std::uint16_t exponent : 5;
        std::uint16_t sign : 1;
    };

    FP16() = default;
    constexpr FP16(std::uint16_t u) : u(u) {};

    FP16(const FP32& fp) {
        *this = fp;
    }

    FP16& operator = (std::uint16_t u){
        this->u = u;
        return *this;
    }

    FP16& operator = (const FP32& fp);

    friend std::ostream& operator<<(std::ostream& os, const FP16& fp);
};

union FP32
{
    std::uint32_t u;
    float f;
    struct
    {
        std::uint32_t mantissa : 23;
        std::uint32_t exponent : 8;
        std::uint32_t sign : 1;
    };

    FP32() = default;
    constexpr FP32(std::uint32_t u) : u(u) {};
    constexpr FP32(float f) : f(f) {};

    FP32(const FP16& fp) {
        *this = fp;
    }

    FP32& operator = (std::uint32_t u){
        this->u = u;
        return *this;
    }

    FP32& operator = (float f) {
        this->f = f;
        return *this;
    }

    FP32& operator = (const FP16& fp);

    friend std::ostream& operator<<(std::ostream& os, const FP32& fp);
};

}