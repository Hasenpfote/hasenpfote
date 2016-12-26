#include <limits>
#include "fp_conversion.h"

namespace hasenpfote{

// Half precision floating point number.
union Half
{
    using binary_type = std::uint16_t;

    static constexpr int MANT_DIG = 11;
    static constexpr int MAX_EXP = 16;
    static constexpr int MIN_EXP = -13;

    static constexpr int BIT_COUNT = 8 * sizeof(binary_type);
    static constexpr int MANTISSA_BIT_COUNT = MANT_DIG - 1;
    static constexpr int EXPONENT_BIT_COUNT = BIT_COUNT - 1 - MANTISSA_BIT_COUNT;

    static constexpr int SIGN_BIT_SHIFT_AMOUNT = BIT_COUNT - 1;
    static constexpr int EXPONENT_BIT_SHIFT_AMOUNT = SIGN_BIT_SHIFT_AMOUNT - EXPONENT_BIT_COUNT;

    static constexpr binary_type SIGN_BIT_MASK = static_cast<binary_type>(1) << SIGN_BIT_SHIFT_AMOUNT;
    static constexpr binary_type MANTISSA_BIT_MASK = static_cast<binary_type>(~static_cast<binary_type>(0)) >> (EXPONENT_BIT_COUNT + 1);
    static constexpr binary_type EXPONENT_BIT_MASK = static_cast<binary_type>(~(SIGN_BIT_MASK | MANTISSA_BIT_MASK));

    static constexpr int EXPONENT_EXCESS = static_cast<binary_type>(~static_cast<binary_type>(0)) >> (BIT_COUNT - EXPONENT_BIT_COUNT + 1);
    static constexpr int BIASED_EXPONENT_UPPER_BOUND = EXPONENT_EXCESS + (MAX_EXP - 1);
    static constexpr int BIASED_EXPONENT_LOWER_BOUND = EXPONENT_EXCESS + (MIN_EXP - 1);

    constexpr explicit Half(binary_type binary) : binary_(binary) {}

    binary_type binary_;
    struct _Bits
    {
        binary_type mantissa_ : MANTISSA_BIT_COUNT;
        binary_type exponent_ : EXPONENT_BIT_COUNT;
        binary_type sign_ : 1;
    } Bits;
};

// Single precision floating point number.
union Single
{
    using raw_type = float;
    using binary_type = std::uint32_t;

    static constexpr int BIT_COUNT = 8 * sizeof(raw_type);
    static constexpr int MANTISSA_BIT_COUNT = std::numeric_limits<raw_type>::digits - 1;
    static constexpr int EXPONENT_BIT_COUNT = BIT_COUNT - 1 - MANTISSA_BIT_COUNT;

    static constexpr int SIGN_BIT_SHIFT_AMOUNT = BIT_COUNT - 1;
    static constexpr int EXPONENT_BIT_SHIFT_AMOUNT = SIGN_BIT_SHIFT_AMOUNT - EXPONENT_BIT_COUNT;

    static constexpr binary_type SIGN_BIT_MASK = static_cast<binary_type>(1) << SIGN_BIT_SHIFT_AMOUNT;
    static constexpr binary_type MANTISSA_BIT_MASK = ~static_cast<binary_type>(0) >> (EXPONENT_BIT_COUNT + 1);
    static constexpr binary_type EXPONENT_BIT_MASK = ~(SIGN_BIT_MASK | MANTISSA_BIT_MASK);
    static constexpr binary_type MANTISSA_HIDDEN_BIT = static_cast<binary_type>(1) << MANTISSA_BIT_COUNT;

    static constexpr int EXPONENT_EXCESS = ~static_cast<binary_type>(0) >> (BIT_COUNT - EXPONENT_BIT_COUNT + 1);

    static constexpr int BIASED_EXPONENT_UPPER_BOUND = EXPONENT_EXCESS + (std::numeric_limits<raw_type>::max_exponent - 1);
    static constexpr int BIASED_EXPONENT_LOWER_BOUND = EXPONENT_EXCESS + (std::numeric_limits<raw_type>::min_exponent - 1);

    constexpr explicit Single(binary_type binary) : binary_(binary) {}
    constexpr explicit Single(raw_type value) : value_(value) {}

    binary_type binary_;
    struct _Bits
    {
        binary_type mantissa_ : MANTISSA_BIT_COUNT;
        binary_type exponent_ : EXPONENT_BIT_COUNT;
        binary_type sign_ : 1;
    } Bits;
    raw_type value_;
};

std::uint16_t ConvertSingleToHalf(float single)
{
    const Single s(single);
    Half h(0);

    /*
        exponent mapping.
            unbiased        -127, -126, ..., -25, ..., -15, -14, ...,   0, ...,  15,  16, ..., 127, 128
            biased(fp32)       0,    1, ..., 102, ..., 112, 113, ..., 127, ..., 142, 143, ..., 254, 255
            biased(fp16)                                 0,   1, ...,  15, ...,  30,  31
    */
    if(s.Bits.exponent_ < Single::BIASED_EXPONENT_LOWER_BOUND){         // exp < 1
        // signed zero or denormal
        h.Bits.exponent_ = 0;   // => signed zero
    }
    else if(s.Bits.exponent_ > Single::BIASED_EXPONENT_UPPER_BOUND){    // exp > 254
        // signed Inf or Nan(SNaN, QNaN)
        h.Bits.exponent_ = Half::BIASED_EXPONENT_UPPER_BOUND + 1;
        h.Bits.mantissa_ = (s.Bits.mantissa_)? 0x200 : 0; // signed NaN => signed QNaN / signed Inf => sigend Inf
    }
    else{                                                               // 1 <= exp <= 254
        // normalized number
        constexpr int REL_EXPONENT_EXCESS = Single::EXPONENT_EXCESS - Half::EXPONENT_EXCESS;
        constexpr int REL_BIASED_EXPONENT_UPPER_BOUND = REL_EXPONENT_EXCESS + Half::BIASED_EXPONENT_UPPER_BOUND;
        constexpr int REL_BIASED_EXPONENT_LOWER_BOUND = REL_EXPONENT_EXCESS + Half::BIASED_EXPONENT_LOWER_BOUND;
        constexpr int REL_BIASED_EXPONENT_DENORM_BOUND = REL_EXPONENT_EXCESS - Half::MANTISSA_BIT_COUNT;

        if(s.Bits.exponent_ > REL_BIASED_EXPONENT_UPPER_BOUND){         // exp > 142
            // overflow
            h.Bits.exponent_ = Half::BIASED_EXPONENT_UPPER_BOUND + 1;   // => signed Inf
        }
        else if((s.Bits.exponent_ < REL_BIASED_EXPONENT_LOWER_BOUND)    // 102 <= exp < 113
              &&(s.Bits.exponent_ >= REL_BIASED_EXPONENT_DENORM_BOUND)){
            // underflow & denormal
            constexpr int REL_MANTISSA_SHIFT_AMOUNT = (Single::MANTISSA_BIT_COUNT + 1) - Half::MANTISSA_BIT_COUNT;
            const auto mantissa = s.Bits.mantissa_ | Single::MANTISSA_HIDDEN_BIT;   // add hidden 1 bit
            const int shift_amount = REL_MANTISSA_SHIFT_AMOUNT + (REL_EXPONENT_EXCESS - s.Bits.exponent_);
            h.Bits.mantissa_ = (mantissa >> shift_amount) & Half::MANTISSA_BIT_MASK;
            if((mantissa >> (shift_amount - 1)) & 0x1)  // round
                h.binary_++;    // 指数ビットにオーバーフローするが問題はない
        }
        else{                                                           // 113 <= exp <= 142
            constexpr int REL_MANTISSA_SHIFT_AMOUNT = Single::MANTISSA_BIT_COUNT - Half::MANTISSA_BIT_COUNT;
            h.Bits.exponent_ = (s.Bits.exponent_ - REL_EXPONENT_EXCESS) & (Half::EXPONENT_BIT_MASK >> Half::EXPONENT_BIT_SHIFT_AMOUNT);
            h.Bits.mantissa_ = (s.Bits.mantissa_ >> REL_MANTISSA_SHIFT_AMOUNT) & Half::MANTISSA_BIT_MASK;
            constexpr auto ROUNDIG_BIT = static_cast<Single::binary_type>(1) << (REL_MANTISSA_SHIFT_AMOUNT - 1);    // 0x1000
            if(s.Bits.mantissa_ & ROUNDIG_BIT)  // round
                h.binary_++;    // 指数ビットにオーバーフローするが問題はない
        }
    }

    h.Bits.sign_ = s.Bits.sign_;

    return h.binary_;
}

float ConvertHalfToSingle(std::uint16_t half)
{
    const Half h(half);
    Single s(static_cast<Single::binary_type>(0));

    constexpr int REL_EXPONENT_EXCESS = Single::EXPONENT_EXCESS - Half::EXPONENT_EXCESS;
    constexpr int REL_MANTISSA_SHIFT_AMOUNT = Single::MANTISSA_BIT_COUNT - Half::MANTISSA_BIT_COUNT;

    s.Bits.mantissa_ = static_cast<Single::binary_type>(h.Bits.mantissa_ << REL_MANTISSA_SHIFT_AMOUNT) & Single::MANTISSA_BIT_MASK;

    if(h.Bits.exponent_ < Half::BIASED_EXPONENT_LOWER_BOUND){
        // signed zero or denormal
        s.Bits.exponent_ = REL_EXPONENT_EXCESS + 1;
        constexpr Single magic(static_cast<Single::binary_type>((REL_EXPONENT_EXCESS + 1) << Single::MANTISSA_BIT_COUNT));
        s.value_ -= magic.value_;   // re-normalize
    }
    else if(h.Bits.exponent_ > Half::BIASED_EXPONENT_UPPER_BOUND){
        // signed Inf or NaN
        s.Bits.exponent_ = Single::BIASED_EXPONENT_UPPER_BOUND + 1;
    }
    else{
        s.Bits.exponent_ = (h.Bits.exponent_ + REL_EXPONENT_EXCESS) & (Single::EXPONENT_BIT_MASK >> Single::EXPONENT_BIT_SHIFT_AMOUNT);
    }

    s.Bits.sign_ = h.Bits.sign_;

    return s.value_;
}

}