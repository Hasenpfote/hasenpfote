#include <cassert>
#include <string>
#include <iostream>
#include <algorithm>
#include "utility.h"
#include "vector2.h"

namespace hasenpfote{ namespace math{

const Vector2 Vector2::ZERO = Vector2(0.0f, 0.0f);
const Vector2 Vector2::E_X = Vector2(1.0f, 0.0f);
const Vector2 Vector2::E_Y = Vector2(0.0f, 1.0f);

Vector2::Vector2(const Vector2& v)
{
    x = v.x;
    y = v.y;
}

Vector2::Vector2(float x, float y)
{
    this->x = x;
    this->y = y;
}

Vector2::Vector2(const std::array<float, 2>& v)
{
    x = v[0];
    y = v[1];
}

Vector2& Vector2::operator = (const Vector2& v)
{
    x = v.x;
    y = v.y;
    return *this;
}

Vector2& Vector2::operator = (const std::array<float, 2>& v)
{
    x = v[0];
    y = v[1];
    return *this;
}

Vector2& Vector2::operator += (const Vector2& v)
{
    x += v.x;
    y += v.y;
    return *this;
}

Vector2& Vector2::operator -= (const Vector2& v)
{
    x -= v.x;
    y -= v.y;
    return *this;
}

Vector2& Vector2::operator *= (float scale)
{
    x *= scale;
    y *= scale;
    return *this;
}

Vector2& Vector2::operator /= (float divisor)
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    x /= divisor;
    y /= divisor;
    return *this;
}

const Vector2 Vector2::operator + () const
{
    return *this;
}

const Vector2 Vector2::operator - () const
{
    return Vector2(-x, -y);
}

const Vector2 Vector2::operator + (const Vector2& v) const
{
    return Vector2(x + v.x, y + v.y);
}

const Vector2 Vector2::operator - (const Vector2& v) const
{
    return Vector2(x - v.x, y - v.y);
}

const Vector2 Vector2::operator * (float scale) const
{
    return Vector2(x * scale, y * scale);
}

const Vector2 Vector2::operator / (float divisor) const
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    return Vector2(x / divisor, y / divisor);
}

const Vector2 operator * (float scale, const Vector2& v)
{
    return Vector2(scale * v.x, scale * v.y);
}

float Vector2::Magnitude() const
{
    return std::sqrtf(MagnitudeSquared());
}

float Vector2::MagnitudeSquared() const
{
    return x * x + y * y;
}

void Vector2::Normalize()
{
    const float mag = Magnitude();
    assert(mag > 0.0f);    // division by zero.
    *this /= mag;
}

Vector2 Vector2::Normalized() const
{
    Vector2 result(*this);
    result.Normalize();
    return result;
}

void Vector2::Negate()
{
    x = -x;
    y = -y;
}

float Vector2::DotProduct(const Vector2& a, const Vector2& b)
{
    return a.x * b.x + a.y * b.y;
}

float Vector2::CrossProduct(const Vector2& a, const Vector2& b)
{
    return a.x * b.y - a.y * b.x;
}

float Vector2::Angle(const Vector2& a, const Vector2& b)
{
    return std::acosf(DotProduct(a, b));
}

Vector2 Vector2::Lerp(const Vector2& a, const Vector2& b, float t)
{
    Vector2 result;
    result.x = a.x + (b.x - a.x) * t;
    result.y = a.y + (b.y - a.y) * t;
    return result;
}

Vector2 Vector2::Slerp(const Vector2& a, const Vector2& b, float t)
{
    Vector2 result;
    const float theta = Angle(a, b);
    if(theta > 0.0f){
        const float fx = std::sinf(theta * (1.0f - t));
        const float fy = std::sinf(theta * t);
        const float cosec = 1.0f / std::sinf(theta);
        result.x = (fx * a.x + fy * b.x) * cosec;
        result.y = (fx * a.y + fy * b.y) * cosec;
    }
    else{
        result = a;
    }
    return result;
}

Vector2 Vector2::BaryCentric(const Vector2& a, const Vector2& b, const Vector2& c, float f, float g)
{
    Vector2 result;
    result.x = a.x + f * (b.x - a.x) + g * (c.x - a.x);
    result.y = a.y + f * (b.y - a.y) + g * (c.y - a.y);
    return result;
}

bool Vector2::IsPerpendicular(const Vector2& a, const Vector2& b)
{
    return !(std::abs(DotProduct(a, b)) > 0.0f);
}

bool Vector2::IsParallel(const Vector2& a, const Vector2& b)
{
    assert(almost_equals(1.0f, a.MagnitudeSquared(), 1));    // a is not an unit vector.
    assert(almost_equals(1.0f, b.MagnitudeSquared(), 1));    // b is not an unit vector.
    return !(std::abs(DotProduct(a, b)) < 1.0f);

}

std::string Vector2::ToString() const
{
    return "Vector2{" + std::to_string(x) + ", " + std::to_string(y) + "}";
}

}}