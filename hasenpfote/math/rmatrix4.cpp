﻿#include <cassert>
#include <string>
#include "utility.h"
#include "vector3.h"
#include "quaternion.h"
#include "rmatrix4.h"

namespace hasenpfote{ namespace math{

const RMatrix4 RMatrix4::IDENTITY = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
};

RMatrix4::RMatrix4(const RMatrix4& m)
{
    *this = m;
}

RMatrix4::RMatrix4(
    float m11, float m12, float m13, float m14,
    float m21, float m22, float m23, float m24,
    float m31, float m32, float m33, float m34,
    float m41, float m42, float m43, float m44)
{
    this->m11 = m11; this->m12 = m12; this->m13 = m13; this->m14 = m14;
    this->m21 = m21; this->m22 = m22; this->m23 = m23; this->m24 = m24;
    this->m31 = m31; this->m32 = m32; this->m33 = m33; this->m34 = m34;
    this->m41 = m41; this->m42 = m42; this->m43 = m43; this->m44 = m44;
}

RMatrix4::RMatrix4(const std::array<float, num_elements>& m)
{
    *this = m;
}

RMatrix4& RMatrix4::operator = (const RMatrix4& m)
{
    std::memcpy(this->m, m.m, sizeof(float) * num_elements);
    return *this;
}

RMatrix4& RMatrix4::operator = (const std::array<float, num_elements>& m)
{
    std::memcpy(this->m, m.data(), sizeof(float) * num_elements);
    return *this;
}

RMatrix4& RMatrix4::operator += (const RMatrix4& m)
{
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] += m.m[i][j];
    return *this;
}

RMatrix4& RMatrix4::operator -= (const RMatrix4& m)
{
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] -= m.m[i][j];
    return *this;
}

RMatrix4& RMatrix4::operator *= (const RMatrix4& m)
{
    *this = *this * m;
    return *this;
}

RMatrix4& RMatrix4::operator *= (float scale)
{
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] *= scale;
    return *this;
}

RMatrix4& RMatrix4::operator /= (float divisor)
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] /= divisor;
    return *this;
}

const RMatrix4 RMatrix4::operator + () const
{
    return *this;
}

const RMatrix4 RMatrix4::operator - () const
{
    RMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = -m[i][j];
    return result;
}

const RMatrix4 RMatrix4::operator + (const RMatrix4& m) const
{
    RMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = this->m[i][j] + m.m[i][j];
    return result;
}

const RMatrix4 RMatrix4::operator - (const RMatrix4& m) const
{
    RMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = this->m[i][j] - m.m[i][j];
    return result;
}

const RMatrix4 RMatrix4::operator * (const RMatrix4& m) const
{
    RMatrix4 result;
    for(auto row = 0; row < order; row++){
        for(auto col = 0; col < order; col++){
            result.m[row][col] = 0.0f;
            for(auto i = 0; i < order; i++){
                result.m[row][col] += (this->m)[row][i] * m.m[i][col];
            }
        }
    }
    return result;
}

const RMatrix4 RMatrix4::operator * (float scale) const
{
    RMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = m[i][j] * scale;
    return result;
}

const RMatrix4 RMatrix4::operator / (float divisor) const
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    RMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = m[i][j] / divisor;
    return result;
}

const RMatrix4 operator * (float scale, const RMatrix4& m)
{
    RMatrix4 result;
    for(auto i = 0; i < RMatrix4::order; i++)
        for(auto j = 0; j < RMatrix4::order; j++)
            result.m[i][j] = scale * m.m[i][j];
    return result;
}

float& RMatrix4::operator () (std::size_t row, std::size_t column)
{
    assert(row < order);
    assert(column < order);
    return m[row][column];
}

const float& RMatrix4::operator () (std::size_t row, std::size_t column) const
{
    assert(row < order);
    assert(column < order);
    return m[row][column];
}

float RMatrix4::Trace() const
{
    return m11 + m22 + m33 + m44;
}

Quaternion RMatrix4::ToRotationQuaternion() const
{
    Quaternion result;
    const float tr = Trace();
    if(tr >= 1.0f){
        // |w| が最大
        result.w = std::sqrtf(tr) * 0.5f;
        const float rcp_4w = 1.0f / (4.0f * result.w);  // 1/4|w|
        result.x = (m23 - m32) * rcp_4w;    // 4wx / 4|w|
        result.y = (m31 - m13) * rcp_4w;    // 4wy / 4|w|
        result.z = (m12 - m21) * rcp_4w;    // 4wz / 4|w|
    }
    else
    if((m11 > m22) && (m11 > m33)){
        // |x| が最大
        result.x = std::sqrtf(m11 - m22 - m33 + 1.0f) * 0.5f;
        const float rcp_4x = 1.0f / (4.0f * result.x);  // 1/4|x|
        result.y = (m21 + m12) * rcp_4x;    // 4xy / 4|x|
        result.z = (m31 + m13) * rcp_4x;    // 4xz / 4|x|
        result.w = (m23 - m32) * rcp_4x;    // 4wx / 4|x|
    }
    else
    if((m22 > m33)){
        // |y| が最大
        result.y = std::sqrtf(m22 - m33 - m11 + 1.0f) * 0.5f;
        const float rcp_4y = 1.0f / (4.0f * result.y);  // 1/4|y|
        result.x = (m21 + m12) * rcp_4y;    // 4xy / 4|y|
        result.z = (m23 + m32) * rcp_4y;    // 4yz / 4|y|
        result.w = (m31 - m13) * rcp_4y;    // 4wy / 4|y|
    }
    else{
        // |z| が最大
        result.z = std::sqrtf(m33 - m11 - m22 + 1.0f) * 0.5f;
        const float rcp_4z = 1.0f / (4.0f * result.z);  // 1/4|z|
        result.x = (m31 + m13) * rcp_4z;    // 4xz / 4|z|
        result.y = (m23 + m32) * rcp_4z;    // 4yz / 4|z|
        result.w = (m12 - m21) * rcp_4z;    // 4wz / 4|z|
    }
    return result;
}

RMatrix4 RMatrix4::Transpose(const RMatrix4& m)
{
    RMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[j][i] = m.m[i][j];
    return result;
}

RMatrix4 RMatrix4::Inverse(const RMatrix4& m, float* determinant)
{
    const float det = m.m11 * (m.m22 * m.m33 - m.m32 * m.m23)
                    - m.m21 * (m.m12 * m.m33 - m.m32 * m.m13)
                    + m.m31 * (m.m12 * m.m23 - m.m22 * m.m13);

    if(determinant)
        *determinant = det;

    if(almost_equals(0.0f, std::fabsf(det), 1))
        return RMatrix4::IDENTITY;

    RMatrix4 result;

    result.m11 = (m.m22 * m.m33 - m.m32 * m.m23) / det;
    result.m12 =-(m.m12 * m.m33 - m.m32 * m.m13) / det;
    result.m13 = (m.m12 * m.m23 - m.m22 * m.m13) / det;
    result.m14 = 0.0f;

    result.m21 =-(m.m21 * m.m33 - m.m31 * m.m23) / det;
    result.m22 = (m.m11 * m.m33 - m.m31 * m.m13) / det;
    result.m23 =-(m.m11 * m.m23 - m.m21 * m.m13) / det;
    result.m24 = 0.0f;

    result.m31 = (m.m21 * m.m32 - m.m31 * m.m22) / det;
    result.m32 =-(m.m11 * m.m32 - m.m31 * m.m12) / det;
    result.m33 = (m.m11 * m.m22 - m.m21 * m.m12) / det;
    result.m34 = 0.0f;

    result.m41 = -(m.m41 * result.m11 + m.m42 * result.m21 + m.m43 * result.m31);
    result.m42 = -(m.m41 * result.m12 + m.m42 * result.m22 + m.m43 * result.m32);
    result.m43 = -(m.m41 * result.m13 + m.m42 * result.m23 + m.m43 * result.m33);
    result.m44 = 1.0f;

    return result;
}

RMatrix4 RMatrix4::Translation(float x, float y, float z)
{
    RMatrix4 result = RMatrix4::IDENTITY;
    result.m41 = x;
    result.m42 = y;
    result.m43 = z;
    return result;
}

RMatrix4 RMatrix4::Scaling(float x, float y, float z)
{
    RMatrix4 result = RMatrix4::IDENTITY;
    result.m11 = x;
    result.m22 = y;
    result.m33 = z;
    return result;
}

RMatrix4 RMatrix4::RotationX(float angle)
{
    RMatrix4 result = RMatrix4::IDENTITY;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    result.m22 = c;
    result.m23 = s;
    result.m32 =-s;
    result.m33 = c;
    return result;
}

RMatrix4 RMatrix4::RotationY(float angle)
{
    RMatrix4 result = RMatrix4::IDENTITY;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    result.m11 = c;
    result.m13 =-s;
    result.m31 = s;
    result.m33 = c;
    return result;
}

RMatrix4 RMatrix4::RotationZ(float angle)
{
    RMatrix4 result = RMatrix4::IDENTITY;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    result.m11 = c;
    result.m12 = s;
    result.m21 =-s;
    result.m22 = c;
    return result;
}

RMatrix4 RMatrix4::RotationAxis(Vector3 axis, float angle)
{
    assert(almost_equals(1.0f, axis.Magnitude(), 1));    // axis is not an unit vector.
    const float x = axis.x;
    const float y = axis.y;
    const float z = axis.z;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    const float vers = 1.0f - c;
    return RMatrix4(
        x*x*vers+c,   x*y*vers+z*s,  x*z*vers-y*s, 0.0f,
        x*y*vers-z*s, y*y*vers+c,    y*z*vers+x*s, 0.0f,
        x*z*vers+y*s, y*z*vers-x*s,  z*z*vers+c,   0.0f,
        0.0f,         0.0f,          0.0f,         1.0f);
}

RMatrix4 RMatrix4::LookAt(Vector3 position, Vector3 target, Vector3 up)
{
    Vector3 zaxis = target - position;
    zaxis.Normalize();
    Vector3 xaxis = Vector3::CrossProduct(up, zaxis);
    xaxis.Normalize();
    Vector3 yaxis = Vector3::CrossProduct(zaxis, xaxis);
    return RMatrix4(
        xaxis.x, yaxis.x, zaxis.x, 0.0f,
        xaxis.y, yaxis.y, zaxis.y, 0.0f,
        xaxis.z, yaxis.z, zaxis.z, 0.0f,
        -Vector3::DotProduct(xaxis, position), -Vector3::DotProduct(yaxis, position), -Vector3::DotProduct(zaxis, position), 1.0f);
}

RMatrix4 RMatrix4::Perspective(float fovy, float aspectRatio, float near, float far)
{
    assert(fovy > 0.0f);
    assert(aspectRatio > 0.0f);
    assert(far > near);
    const float cot = 1.0f / std::tanf(fovy * 0.5f);
    const float q = far / (far - near);
    return RMatrix4(
        cot/aspectRatio, 0.0f, 0.0f,    0.0f,
        0.0f,            cot,  0.0f,    0.0f,
        0.0f,            0.0f, q,       1.0f,
        0.0f,            0.0f, -near*q, 0.0f);
}

RMatrix4 RMatrix4::Frustum(float top, float bottom, float left, float right, float near, float far)
{
    assert(top > bottom);
    assert(right > left);
    assert(far > near);
    const float w = 2.0f * near / (right - left);
    const float h = 2.0f * near / (top - bottom);
    const float q = far / (far - near);
    const float woff = -(right + left) / (right - left);
    const float hoff = -(top + bottom) / (top - bottom);
    return RMatrix4(
        w,    0.0f, 0.0f,    0.0f,
        0.0f, h,    0.0f,    0.0f,
        woff, hoff, q,       1.0f,
        0.0f, 0.0f, -near*q, 0.0f);
}

RMatrix4 RMatrix4::Ortho(float top, float bottom, float left, float right, float near, float far)
{
    assert(top > bottom);
    assert(right > left);
    assert(far > near);
    const float w = 2.0f / (right - left);
    const float h = 2.0f / (top - bottom);
    const float q = 1.0f / (far - near);
    const float woff = -(right + left) / (right - left);
    const float hoff = -(top + bottom) / (top - bottom);
    return RMatrix4(
        w,    0.0f, 0.0f,    0.0f,
        0.0f, h,    0.0f,    0.0f,
        0.0f, 0.0f, q,       0.0f,
        woff, hoff, -near*q, 1.0f);
}

std::ostream& operator<<(std::ostream& os, const RMatrix4& m)
{
    const auto flags = os.flags();
    for(auto i = 0; i < RMatrix4::order; i++){
        os << "RMatrix4[" << i << "]{" << m.m[i][0];
        for (auto j = 1; j < RMatrix4::order; j++){
            os << ", " << m.m[i][j];
        }
        os << "}" << std::endl;
    }
    os.flags(flags);
    return os;
}

}}