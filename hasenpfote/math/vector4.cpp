#include "../assert.h"
#include <string>
#include <iostream>
#include <algorithm>
#include "utility.h"
#include "vector3.h"
#include "cmatrix4.h"
#include "rmatrix4.h"
#include "vector4.h"

namespace hasenpfote{ namespace math{

const Vector4 Vector4::ZERO = Vector4(0.0f, 0.0f, 0.0f, 0.0f);

Vector4::Vector4(const Vector4& v)
    : x(v.x), y(v.y), z(v.z), w(v.w)
{
}

Vector4::Vector4(float x, float y, float z, float w)
    : x(x), y(y), z(z), w(w)
{
}

Vector4::Vector4(const Vector3& v, float w)
    : Vector4(v.GetX(), v.GetY(), v.GetZ(), w)
{
   
}

Vector4::Vector4(const Array& v)
    : Vector4(v[0], v[1], v[2], v[3])
{
}

Vector4& Vector4::operator = (const Vector4& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = v.w;
    return *this;
}

Vector4& Vector4::operator = (const Array& v)
{
    x = v[0];
    y = v[1];
    z = v[2];
    w = v[3];
    return *this;
}

Vector4& Vector4::operator += (const Vector4& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;
    return *this;
}

Vector4& Vector4::operator -= (const Vector4& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;
    return *this;
}

Vector4& Vector4::operator *= (float scale)
{
    x *= scale;
    y *= scale;
    z *= scale;
    w *= scale;
    return *this;
}

Vector4& Vector4::operator /= (float divisor)
{
    ASSERT_MSG(std::fabsf(divisor) > 0.0f, "Division by zero.");
    x /= divisor;
    y /= divisor;
    z /= divisor;
    w /= divisor;
    return *this;
}

float Vector4::Magnitude() const
{
    return std::sqrtf(MagnitudeSquared());
}

float Vector4::MagnitudeSquared() const
{
    return x * x + y * y + z * z + w * w;
}

void Vector4::Normalize()
{
    const float mag = Magnitude();
    ASSERT_MSG(mag > 0.0f, "Division by zero.");
    *this /= mag;
}

Vector4 Vector4::Normalized() const
{
    Vector4 result(*this);
    result.Normalize();
    return result;
}

void Vector4::Negate()
{
    x = -x;
    y = -y;
    z = -z;
    w = -w;
}

float Vector4::DotProduct(const Vector4& a, const Vector4& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

Vector4 Vector4::Minimize(const Vector4& a, const Vector4& b)
{
    return Vector4(
        (a.x < b.x)? a.x : b.x,
        (a.y < b.y)? a.y : b.y,
        (a.z < b.z)? a.z : b.z,
        (a.w < b.w)? a.w : b.w
    );
}

Vector4 Vector4::Maximize(const Vector4& a, const Vector4& b)
{
    return Vector4(
        (a.x > b.x)? a.x : b.x,
        (a.y > b.y)? a.y : b.y,
        (a.z > b.z)? a.z : b.z,
        (a.w > b.w)? a.w : b.w
    );
}

Vector4 operator + (const Vector4& v)
{
    return v;
}

Vector4 operator - (const Vector4& v)
{
    return Vector4(-v.GetX(), -v.GetY(), -v.GetZ(), -v.GetW());
}

Vector4 operator + (const Vector4& lhs, const Vector4& rhs)
{
    return Vector4(lhs) += rhs;
}

Vector4 operator - (const Vector4& lhs, const Vector4& rhs)
{
    return Vector4(lhs) -= rhs;
}

Vector4 operator * (const Vector4& v, float scale)
{
    return Vector4(v) *= scale;
}

Vector4 operator * (float scale, const Vector4& v)
{
    return Vector4(v) *= scale;
}

Vector4 operator / (const Vector4& v, float divisor)
{
    return Vector4(v) /= divisor;
}

Vector4 operator * (const CMatrix4& m, const Vector4& v)
{
    return Vector4(
        Vector4::DotProduct(m.GetRow(0), v),
        Vector4::DotProduct(m.GetRow(1), v),
        Vector4::DotProduct(m.GetRow(2), v),
        Vector4::DotProduct(m.GetRow(3), v)
    );
}

Vector4 operator * (const Vector4& v, const RMatrix4& m)
{
    return Vector4(
        Vector4::DotProduct(v, m.GetColumn(0)),
        Vector4::DotProduct(v, m.GetColumn(1)),
        Vector4::DotProduct(v, m.GetColumn(2)),
        Vector4::DotProduct(v, m.GetColumn(3))
    );
}

std::ostream& operator << (std::ostream& os, const Vector4& v)
{
    const auto flags = os.flags();
    os << "Vector4{" << v.GetX() << ", " << v.GetY() << ", " << v.GetZ() << ", " << v.GetW() << "}";
    os.flags(flags);
    return os;
}

}}