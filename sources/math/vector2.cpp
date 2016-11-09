#include <cstdint>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include "../assert.h"
#include "utility.h"
#include "vector2.h"

namespace hasenpfote{ namespace math{

const Vector2 Vector2::ZERO = Vector2(0.0f, 0.0f);
const Vector2 Vector2::E_X = Vector2(1.0f, 0.0f);
const Vector2 Vector2::E_Y = Vector2(0.0f, 1.0f);

Vector2::Vector2(const Vector2& v)
    : x(v.x), y(v.y)
{
}

Vector2::Vector2(float x, float y)
    : x(x), y(y)
{
}

Vector2::Vector2(const Array& v)
    : Vector2(v[0], v[1])
{
}

Vector2& Vector2::operator = (const Vector2& v)
{
    x = v.x;
    y = v.y;
    return *this;
}

Vector2& Vector2::operator = (const Array& v)
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
    HASENPFOTE_ASSERT_MSG(std::abs(divisor) > 0.0f, "Division by zero.");
    x /= divisor;
    y /= divisor;
    return *this;
}

float Vector2::Magnitude() const
{
    return std::sqrt(MagnitudeSquared());
}

float Vector2::MagnitudeSquared() const
{
    return x * x + y * y;
}

void Vector2::Normalize()
{
    const float mag = Magnitude();
    HASENPFOTE_ASSERT_MSG(mag > 0.0f, "Division by zero.");
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
    return std::acos(DotProduct(a, b));
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
        const float fx = std::sin(theta * (1.0f - t));
        const float fy = std::sin(theta * t);
        const float cosec = 1.0f / std::sin(theta);
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
    HASENPFOTE_ASSERT_MSG(almost_equals(1.0f, a.MagnitudeSquared(), 1), "Not an unit vector.");
    HASENPFOTE_ASSERT_MSG(almost_equals(1.0f, b.MagnitudeSquared(), 1), "Not an unit vector.");
    return !(std::abs(DotProduct(a, b)) < 1.0f);
}

Vector2 Vector2::Minimize(const Vector2& a, const Vector2& b)
{
    return Vector2(
        (a.x < b.x)? a.x : b.x,
        (a.y < b.y)? a.y : b.y
    );
}

Vector2 Vector2::Maximize(const Vector2& a, const Vector2& b)
{
    return Vector2(
        (a.x > b.x)? a.x : b.x,
        (a.y > b.y)? a.y : b.y
    );
}

Vector2 operator + (const Vector2& v)
{
    return v;
}

Vector2 operator - (const Vector2& v)
{
    return Vector2(-v.GetX(), -v.GetY());
}

Vector2 operator + (const Vector2& lhs, const Vector2& rhs)
{
    return Vector2(lhs) += rhs;
}

Vector2 operator - (const Vector2& lhs, const Vector2& rhs)
{
    return Vector2(lhs) -= rhs;
}

Vector2 operator * (const Vector2& v, float scale)
{
    return Vector2(v) *= scale;
}

Vector2 operator * (float scale, const Vector2& v)
{
    return Vector2(v) *= scale;
}

Vector2 operator / (const Vector2& v, float divisor)
{
    return Vector2(v) /= divisor;
}

std::ostream& operator << (std::ostream& os, const Vector2& v)
{
    const auto flags = os.flags();
    os << "Vector2{" << v.GetX() << ", " << v.GetY() << "}";
    os.flags(flags);
    return os;
}

}}