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
    HASENPFOTE_ASSERT_MSG(std::fabsf(divisor) > 0.0f, "Division by zero.");
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
    return std::sqrtf(NormSquared());
}

float ComplexNumber::Argument() const
{
    return std::atan2f(im, re);
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
    return ComplexNumber(rho * std::cosf(theta), rho * std::sinf(theta));
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
    const float rho = std::powf(c.Norm(), exponent);
    const float arg = c.Argument() * exponent;
    return Polar(rho, arg);
}

ComplexNumber ComplexNumber::Ln(const ComplexNumber& c)
{
    return ComplexNumber(std::logf(c.Norm()), c.Argument());
}

ComplexNumber ComplexNumber::Exp(const ComplexNumber& c)
{
    return Polar(std::expf(c.re), c.im);
}

ComplexNumber ComplexNumber::Rotation(float angle)
{
    return ComplexNumber(std::cosf(angle), std::sinf(angle));
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