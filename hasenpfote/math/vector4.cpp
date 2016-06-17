#include <cassert>
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
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = v.w;
}

Vector4::Vector4(float x, float y, float z, float w)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

Vector4::Vector4(const Vector3& v, float w)
{
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
    this->w = w;
}

Vector4::Vector4(const std::array<float, 4>& v)
{
    x = v[0];
    y = v[1];
    z = v[2];
    w = v[3];
}

Vector4& Vector4::operator = (const Vector4& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = v.w;
    return *this;
}

Vector4& Vector4::operator = (const std::array<float, 4>& v)
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
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    x /= divisor;
    y /= divisor;
    z /= divisor;
    w /= divisor;
    return *this;
}

const Vector4 Vector4::operator + () const
{
    return *this;
}

const Vector4 Vector4::operator - () const
{
    return Vector4(-x, -y, -z, -w);
}

const Vector4 Vector4::operator + (const Vector4& v) const
{
    return Vector4(x + v.x, y + v.y, z + v.z, w + v.w);
}

const Vector4 Vector4::operator - (const Vector4& v) const
{
    return Vector4(x - v.x, y - v.y, z - v.z, w - v.w);
}

const Vector4 Vector4::operator * (float scale) const
{
    return Vector4(x * scale, y * scale, z * scale, w * scale);
}

const Vector4 Vector4::operator / (float divisor) const
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    return Vector4(x / divisor, y / divisor, z / divisor, w / divisor);
}

const Vector4 operator * (float scale, const Vector4& v)
{
    return Vector4(scale * v.x, scale * v.y, scale * v.z, scale * v.w);
}

const Vector4 operator * (const CMatrix4& m, const Vector4& v)
{
    Vector4 result;
    result.x = m.m11 * v.x + m.m12 * v.y + m.m13 * v.z + m.m14 * v.w;
    result.y = m.m21 * v.x + m.m22 * v.y + m.m23 * v.z + m.m24 * v.w;
    result.z = m.m31 * v.x + m.m32 * v.y + m.m33 * v.z + m.m34 * v.w;
    result.w = m.m41 * v.x + m.m42 * v.y + m.m43 * v.z + m.m44 * v.w;
    return result;
}

const Vector4 operator * (const Vector4& v, const RMatrix4& m)
{
    Vector4 result;
    result.x = v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + v.w * m.m41;
    result.y = v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + v.w * m.m42;
    result.z = v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + v.w * m.m43;
    result.w = v.x * m.m14 + v.y * m.m24 + v.z * m.m34 + v.w * m.m44;
    return result;
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
    assert(mag > 0.0f);    // division by zero.
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

std::ostream& operator<<(std::ostream& os, const Vector4& v)
{
    const auto flags = os.flags();
    os << "Vector4{" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "}";
    os.flags(flags);
    return os;
}

}}