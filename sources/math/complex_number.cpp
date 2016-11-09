#include <cmath>
#include "../assert.h"
#include "complex_number.h"

namespace hasenpfote{ namespace math{

const ComplexNumber ComplexNumber::IDENTITY = ComplexNumber(1.0f, 0.0f);

ComplexNumber::ComplexNumber(const ComplexNumber& c)
    : re(c.re), im(c.im)
{
}

ComplexNumber::ComplexNumber(float re, float im)
    : re(re), im(im)
{
}

ComplexNumber& ComplexNumber::operator = (const ComplexNumber& c)
{
    re = c.re;
    im = c.im;
    return *this;
}

ComplexNumber& ComplexNumber::operator += (const ComplexNumber& c)
{
    re += c.re;
    im += c.im;
    return *this;
}

ComplexNumber& ComplexNumber::operator -= (const ComplexNumber& c)
{
    re -= c.re;
    im -= c.im;
    return *this;
}

ComplexNumber& ComplexNumber::operator *= (const ComplexNumber& c)
{
    const ComplexNumber temp(*this);
    re = temp.re * c.re - temp.im * c.im;
    im = temp.re * c.im + temp.im * c.re;
    return *this;
}

ComplexNumber& ComplexNumber::operator *= (float scale)
{
    re *= scale;
    im *= scale;
    return *this;
}

ComplexNumber& ComplexNumber::operator /= (float divisor)
{
    HASENPFOTE_ASSERT_MSG(std::abs(divisor) > 0.0f, "Division by zero.");
    re /= divisor;
    im /= divisor;
    return *this;
}

float ComplexNumber::NormSquared() const
{
    return re * re + im * im;
}

float ComplexNumber::Norm() const
{
    return std::sqrt(NormSquared());
}

float ComplexNumber::Argument() const
{
    return std::atan2(im, re);
}

void ComplexNumber::Normalize()
{
    const float n = Norm();
    HASENPFOTE_ASSERT_MSG(n > 0.0f, "Division by zero.");
    re /= n;
    im /= n;
}

ComplexNumber ComplexNumber::Normalized() const
{
    const float n = Norm();
    HASENPFOTE_ASSERT_MSG(n > 0.0f, "Division by zero.");
    return ComplexNumber(re / n, im / n);
}

ComplexNumber ComplexNumber::Polar(float rho, float theta)
{
    return ComplexNumber(rho * std::cos(theta), rho * std::sin(theta));
}

ComplexNumber ComplexNumber::Inverse(const ComplexNumber& c)
{
    const float nsq = c.NormSquared();
    HASENPFOTE_ASSERT_MSG(nsq > 0.0f, "Division by zero.");
    return ComplexNumber(c.re / nsq, -c.im / nsq);
}

ComplexNumber ComplexNumber::Conjugate(const ComplexNumber& c)
{
    return ComplexNumber(c.re, -c.im);
}

ComplexNumber ComplexNumber::Pow(const ComplexNumber& c, float exponent)
{
    const float rho = std::pow(c.Norm(), exponent);
    const float arg = c.Argument() * exponent;
    return Polar(rho, arg);
}

ComplexNumber ComplexNumber::Ln(const ComplexNumber& c)
{
    return ComplexNumber(std::log(c.Norm()), c.Argument());
}

ComplexNumber ComplexNumber::Exp(const ComplexNumber& c)
{
    return Polar(std::exp(c.re), c.im);
}

ComplexNumber ComplexNumber::Rotation(float angle)
{
    return ComplexNumber(std::cos(angle), std::sin(angle));
}

ComplexNumber operator + (const ComplexNumber& c)
{
    return c;
}

ComplexNumber operator - (const ComplexNumber& c)
{
    return ComplexNumber(-c.GetRealPart(), -c.GetImaginaryPart());
}

ComplexNumber operator + (const ComplexNumber& lhs, const ComplexNumber& rhs)
{
    return ComplexNumber(lhs) += rhs;
}

ComplexNumber operator - (const ComplexNumber& lhs, const ComplexNumber& rhs)
{
    return ComplexNumber(lhs) -= rhs;
}

ComplexNumber operator * (const ComplexNumber& lhs, const ComplexNumber& rhs)
{
    return ComplexNumber(lhs) *= rhs;
}

ComplexNumber operator * (const ComplexNumber& c, float scale)
{
    return ComplexNumber(c) *= scale;
}

ComplexNumber operator * (float scale, const ComplexNumber& c)
{
    return ComplexNumber(c) *= scale;
}

ComplexNumber operator / (const ComplexNumber& c, float divisor)
{
    return ComplexNumber(c) /= divisor;
}

std::ostream& operator << (std::ostream& os, const ComplexNumber& c)
{
    const auto flags = os.flags();
    os << "ComplexNumber{" << c.GetRealPart() << ", " << c.GetImaginaryPart() << "}";
    os.flags(flags);
    return os;
}

}}