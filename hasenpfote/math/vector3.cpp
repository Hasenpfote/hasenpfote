#include <cassert>
#include <string>
#include <iostream>
#include <algorithm>
#include "utility.h"
#include "cmatrix4.h"
#include "rmatrix4.h"
#include "vector3.h"

namespace hasenpfote{ namespace math{

const Vector3 Vector3::ZERO = Vector3(0.0f, 0.0f, 0.0f);
const Vector3 Vector3::E_X = Vector3(1.0f, 0.0f, 0.0f);
const Vector3 Vector3::E_Y = Vector3(0.0f, 1.0f, 0.0f);
const Vector3 Vector3::E_Z = Vector3(0.0f, 0.0f, 1.0f);

Vector3::Vector3(const Vector3& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

Vector3::Vector3(float x, float y, float z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector3::Vector3(const std::array<float, 3>& v)
{
    x = v[0];
    y = v[1];
    z = v[2];
}

Vector3& Vector3::operator = (const Vector3& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
}

Vector3& Vector3::operator = (const std::array<float, 3>& v)
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
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    x /= divisor;
    y /= divisor;
    z /= divisor;
    return *this;
}

const Vector3 Vector3::operator + () const
{
    return *this;
}

const Vector3 Vector3::operator - () const
{
    return Vector3(-x, -y, -z);
}

const Vector3 Vector3::operator + (const Vector3& v) const
{
    return Vector3(x + v.x, y + v.y, z + v.z);
}

const Vector3 Vector3::operator - (const Vector3& v) const
{
    return Vector3(x - v.x, y - v.y, z - v.z);
}

const Vector3 Vector3::operator * (float scale) const
{
    return Vector3(x * scale, y * scale, z * scale);
}

const Vector3 Vector3::operator / (float divisor) const
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    return Vector3(x / divisor, y / divisor, z / divisor);
}

const Vector3 operator * (float scale, const Vector3& v)
{
    return Vector3(scale * v.x, scale * v.y, scale * v.z);
}

const Vector3 operator * (const CMatrix4& m, const Vector3& v)
{
    Vector3 result;
    result.x = m.m11 * v.x + m.m12 * v.y + m.m13 * v.z + m.m14;
    result.y = m.m21 * v.x + m.m22 * v.y + m.m23 * v.z + m.m24;
    result.z = m.m31 * v.x + m.m32 * v.y + m.m33 * v.z + m.m34;
    return result;
}

const Vector3 operator * (const Vector3& v, const RMatrix4& m)
{
    Vector3 result;
    result.x = v.x * m.m11 + v.y * m.m21 + v.z * m.m31 + m.m41;
    result.y = v.x * m.m12 + v.y * m.m22 + v.z * m.m32 + m.m42;
    result.z = v.x * m.m13 + v.y * m.m23 + v.z * m.m33 + m.m43;
    return result;
}

float Vector3::Magnitude() const
{
    return std::sqrtf(MagnitudeSquared());
}

float Vector3::MagnitudeSquared() const
{
    return x * x + y * y + z * z;
}

void Vector3::Normalize()
{
    const float mag = Magnitude();
    assert(mag > 0.0f);    // division by zero.
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
    return std::acosf(DotProduct(a, b));
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
        const float fx = std::sinf(theta * (1.0f - t));
        const float fy = std::sinf(theta * t);
        const float cosec = 1.0f / std::sinf(theta);
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
    assert(almost_equals(1.0f, a.MagnitudeSquared(), 1));    // a is not an unit vector.
    assert(almost_equals(1.0f, b.MagnitudeSquared(), 1));    // b is not an unit vector.
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

std::ostream& operator<<(std::ostream& os, const Vector3& v)
{
    const auto flags = os.flags();
    os << "Vector3{" << v.x << ", " << v.y << ", " << v.z << "}";
    os.flags(flags);
    return os;
}

}}