#include <cassert>
#include "complex_number.h"

namespace hasenpfote{ namespace math{

const ComplexNumber ComplexNumber::IDENTITY = ComplexNumber(1.0f, 0.0f);

ComplexNumber::ComplexNumber(const ComplexNumber& c)
{
    re = c.re;
    im = c.im;
}

ComplexNumber::ComplexNumber(float re, float im)
{
    this->re = re;
    this->im = im;
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
    *this = *this * c;
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
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    re /= divisor;
    im /= divisor;
    return *this;
}

const ComplexNumber ComplexNumber::operator + () const
{
    return *this;
}

const ComplexNumber ComplexNumber::operator - () const
{
    return ComplexNumber(-re, -im);
}

const ComplexNumber ComplexNumber::operator + (const ComplexNumber& c) const
{
    return ComplexNumber(re + c.re, im + c.im);
}

const ComplexNumber ComplexNumber::operator - (const ComplexNumber& c) const
{
    return ComplexNumber(re - c.re, im - c.im);
}

const ComplexNumber ComplexNumber::operator * (const ComplexNumber& c) const
{
    ComplexNumber result;
    result.re = re * c.re - im * c.im;
    result.im = re * c.im + im * c.re;
    return result;
}

const ComplexNumber ComplexNumber::operator * (float scale) const
{
    return ComplexNumber(re * scale, im * scale);
}

const ComplexNumber ComplexNumber::operator / (float divisor) const
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    return ComplexNumber(re / divisor, im / divisor);
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
    assert(n > 0.0f);    // division by zero.
    re /= n;
    im /= n;
}

ComplexNumber ComplexNumber::Normalized() const
{
    const float n = Norm();
    assert(n > 0.0f);    // division by zero.
    return ComplexNumber(re / n, im / n);
}

ComplexNumber ComplexNumber::Polar(float rho, float theta)
{
    return ComplexNumber(rho * std::cosf(theta), rho * std::sinf(theta));
}

ComplexNumber ComplexNumber::Inverse(const ComplexNumber& c)
{
    const float nsq = c.NormSquared();
    assert(nsq > 0.0f);    // division by zero.
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

std::string ComplexNumber::ToString() const
{
    return "ComplexNumber{" + std::to_string(re) + ", " + std::to_string(im) + "}";
}

}}