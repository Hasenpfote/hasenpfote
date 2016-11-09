#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include "../assert.h"
#include "utility.h"
#include "cmatrix4.h"
#include "rmatrix4.h"
#include "vector4.h"
#include "vector3.h"

namespace hasenpfote{ namespace math{

const Vector3 Vector3::ZERO = Vector3(0.0f, 0.0f, 0.0f);
const Vector3 Vector3::E_X = Vector3(1.0f, 0.0f, 0.0f);
const Vector3 Vector3::E_Y = Vector3(0.0f, 1.0f, 0.0f);
const Vector3 Vector3::E_Z = Vector3(0.0f, 0.0f, 1.0f);

Vector3::Vector3(const Vector3& v)
    : x(v.x), y(v.y), z(v.z)
{
}

Vector3::Vector3(float x, float y, float z)
    : x(x), y(y), z(z)
{
}

Vector3::Vector3(const Array& v)
    : Vector3(v[0], v[1], v[2])
{
}

Vector3::Vector3(const Vector4& v)
    : Vector3(v.GetX(), v.GetY(), v.GetZ())
{
}

Vector3& Vector3::operator = (const Vector3& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
}

Vector3& Vector3::operator = (const Array& v)
{
    x = v[0];
    y = v[1];
    z = v[2];
    return *this;
}

Vector3& Vector3::operator += (const Vector3& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

Vector3& Vector3::operator -= (const Vector3& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

Vector3& Vector3::operator *= (float scale)
{
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;
}

Vector3& Vector3::operator /= (float divisor)
{
    HASENPFOTE_ASSERT_MSG(std::abs(divisor) > 0.0f, "Division by zero.");
    x /= divisor;
    y /= divisor;
    z /= divisor;
    return *this;
}

float Vector3::Magnitude() const
{
    return std::sqrt(MagnitudeSquared());
}

float Vector3::MagnitudeSquared() const
{
    return x * x + y * y + z * z;
}

void Vector3::Normalize()
{
    const float mag = Magnitude();
    HASENPFOTE_ASSERT_MSG(mag > 0.0f, "Division by zero.");
    *this /= mag;
}

Vector3 Vector3::Normalized() const
{
    Vector3 result(*this);
    result.Normalize();
    return result;
}

void Vector3::Negate()
{
    x = -x;
    y = -y;
    z = -z;
}

float Vector3::DotProduct(const Vector3& a, const Vector3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 Vector3::CrossProduct(const Vector3& a, const Vector3& b)
{
    Vector3 result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

float Vector3::Angle(const Vector3& a, const Vector3& b)
{
    return std::acos(DotProduct(a, b));
}

Vector3 Vector3::Lerp(const Vector3& a, const Vector3& b, float t)
{
    Vector3 result;
    result.x = a.x + (b.x - a.x) * t;
    result.y = a.y + (b.y - a.y) * t;
    result.z = a.z + (b.z - a.z) * t;
    return result;
}

Vector3 Vector3::Slerp(const Vector3& a, const Vector3& b, float t)
{
    Vector3 result;
    const float theta = Angle(a, b);
    if(theta > 0.0f){
        const float fx = std::sin(theta * (1.0f - t));
        const float fy = std::sin(theta * t);
        const float cosec = 1.0f / std::sin(theta);
        result.x = (fx * a.x + fy * b.x) * cosec;
        result.y = (fx * a.y + fy * b.y) * cosec;
        result.z = (fx * a.z + fy * b.z) * cosec;
    }
    else{
        result = a;
    }
    return result;
}

Vector3 Vector3::BaryCentric(const Vector3& a, const Vector3& b, const Vector3& c, float f, float g)
{
    Vector3 result;
    result.x = a.x + f * (b.x - a.x) + g * (c.x - a.x);
    result.y = a.y + f * (b.y - a.y) + g * (c.y - a.y);
    result.z = a.z + f * (b.z - a.z) + g * (c.z - a.z);
    return result;
}

bool Vector3::IsPerpendicular(const Vector3& a, const Vector3& b)
{
    return !(std::abs(DotProduct(a, b)) > 0.0f);
}

bool Vector3::IsParallel(const Vector3& a, const Vector3& b)
{
    HASENPFOTE_ASSERT_MSG(almost_equals(1.0f, a.MagnitudeSquared(), 1), "Not an unit vector.");
    HASENPFOTE_ASSERT_MSG(almost_equals(1.0f, b.MagnitudeSquared(), 1), "Not an unit vector.");
    return !(std::abs(DotProduct(a, b)) < 1.0f);
}

Vector3 Vector3::Minimize(const Vector3& a, const Vector3& b)
{
    return Vector3(
        (a.x < b.x)? a.x : b.x,
        (a.y < b.y)? a.y : b.y,
        (a.z < b.z)? a.z : b.z
    );
}

Vector3 Vector3::Maximize(const Vector3& a, const Vector3& b)
{
    return Vector3(
        (a.x > b.x)? a.x : b.x,
        (a.y > b.y)? a.y : b.y,
        (a.z > b.z)? a.z : b.z
    );
}

Vector3 operator + (const Vector3& v)
{
    return v;
}

Vector3 operator - (const Vector3& v)
{
    return Vector3(-v.GetX(), -v.GetY(), -v.GetZ());
}

Vector3 operator + (const Vector3& lhs, const Vector3& rhs)
{
    return Vector3(lhs) += rhs;
}

Vector3 operator - (const Vector3& lhs, const Vector3& rhs)
{
    return Vector3(lhs) -= rhs;
}

Vector3 operator * (const Vector3& v, float scale)
{
    return Vector3(v) *= scale;
}

Vector3 operator * (float scale, const Vector3& v)
{
    return Vector3(v) *= scale;
}

Vector3 operator / (const Vector3& v, float divisor)
{
    return Vector3(v) /= divisor;
}

Vector3 operator * (const CMatrix4& m, const Vector3& v)
{
    const Vector4 h(v, 1.0f);
    return Vector3(
        Vector4::DotProduct(m.GetRow(0), h),
        Vector4::DotProduct(m.GetRow(1), h),
        Vector4::DotProduct(m.GetRow(2), h)
    );
}

Vector3 operator * (const Vector3& v, const RMatrix4& m)
{
    const Vector4 h(v, 1.0f);
    return Vector3(
        Vector4::DotProduct(h, m.GetColumn(0)),
        Vector4::DotProduct(h, m.GetColumn(1)),
        Vector4::DotProduct(h, m.GetColumn(2))
    );
}

std::ostream& operator << (std::ostream& os, const Vector3& v)
{
    const auto flags = os.flags();
    os << "Vector3{" << v.GetX() << ", " << v.GetY() << ", " << v.GetZ() << "}";
    os.flags(flags);
    return os;
}

}}