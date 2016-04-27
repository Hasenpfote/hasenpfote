#include <iomanip>
#include "fp.h"

namespace hasenpfote{

FP16& FP16::operator = (const FP32& fp)
{
    FP16 o = { 0 };
    if(fp.exponent == 0){   // signed zero/denormal
        o.exponent = 0;     // => signed zero
    }
    else if(fp.exponent == 255){    // signed Inf/NaN(SNaN or QNaN)
        o.exponent = 31;
        o.mantissa = (fp.mantissa)? 0x200 : 0; // signed NaN => signed QNaN / signed Inf => sigend Inf
    }
    else{
        // normalized number
        const std::int32_t newexp = fp.exponent - 127 + 15;
        if(newexp >= 31){       // overflow
            o.exponent = 31;    // => signed Inf
        }
        else if(newexp <= 0){   // underflow
            if(newexp >= -10){  
                // 非正規化数の処理(2^-25 <= newexp <= 2^0)
                std::uint32_t mant = fp.mantissa | 0x800000; // hidden 1 bit
                o.mantissa = mant >> (14 - newexp);
                if((mant >> (13 - newexp)) & 1) // round
                    o.u++;  // 指数ビットにオーバーフローする可能性があるが問題ない
            }
        }
        else {
            o.exponent = newexp;
            o.mantissa = fp.mantissa >> 13;
            if(fp.mantissa & 0x1000)    // round
                o.u++;  // 指数ビットにオーバーフローする可能性があるが問題ない
        }
    }
    o.sign = fp.sign;

    this->u = o.u;
    return *this;
}

std::ostream& operator<<(std::ostream& os, const FP16& fp)
{
    const auto flags = os.flags();
    os << "FP16{s=" << fp.sign;
    os << std::setfill('0');
    os << " e=0x" << std::hex << std::setw(2) << fp.exponent << "(" << std::dec << fp.exponent << ")";
    os << " m=0x" << std::hex << std::setw(3) << fp.mantissa << "(" << std::dec << fp.mantissa << ")";
    os << " u=0x" << std::hex << std::setw(4) << fp.u << "}";
    os.flags(flags);
    return os;
}

FP32& FP32::operator = (const FP16& fp)
{
    constexpr FP32 magic = { static_cast<uint32_t>(113 << 23) };    // sign:0 exp:113 mantissa:0
    constexpr std::uint32_t shifted_exp = 0x7c00 << 13;
    FP32 o;
    o.u = (fp.u & 0x7fff) << 13;    // exponent / mantissa bits
    const std::uint32_t exp = o.u & shifted_exp;
    o.u += (127 - 15) << 23;    // FP16 と FP32 のバイアスの差分(127 - 15)を足しこみ指数部を調整

    if(exp == shifted_exp){ // signed Inf/NaN
        o.u += (128 - 16) << 23;    // FP16 と FP32 の INF の差分(128 - 16)を足しこみ指数部を調整
    }
    else if(exp == 0){      // signed zero/denormal
        // FP16 の非正規化数を再正規化する
        o.u += 1 << 23;     // extra exp adjust
        o.f -= magic.f;     // renormalize
    }
    o.u |= (fp.u & 0x8000) << 16;   // sign bit

    this->u = o.u;
    return *this;
}

std::ostream& operator<<(std::ostream& os, const FP32& fp)
{
    const auto flags = os.flags();
    os << "FP32{s=" << fp.sign;
    os << std::setfill('0');
    os << " e=0x" << std::hex << std::setw(2) << fp.exponent << "(" << std::dec << fp.exponent << ")";
    os << " m=0x" << std::hex << std::setw(6) << fp.mantissa << "(" << std::dec << fp.mantissa << ")";
    os << " u=0x" << std::hex << std::setw(8) << fp.u << "}";
    os << " f=" << std::fixed << fp.f << "}";
    os.flags(flags);
    return os;
}

}