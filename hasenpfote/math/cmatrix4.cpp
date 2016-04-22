#include <cassert>
#include <string>
#include "utility.h"
#include "vector3.h"
#include "quaternion.h"
#include "cmatrix4.h"

namespace hasenpfote{ namespace math{

const CMatrix4 CMatrix4::IDENTITY = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
};

CMatrix4::CMatrix4(const CMatrix4& m)
{
    *this = m;
}

CMatrix4::CMatrix4(
    float m11, float m12, float m13, float m14,
    float m21, float m22, float m23, float m24,
    float m31, float m32, float m33, float m34,
    float m41, float m42, float m43, float m44)
{
    this->m11 = m11; this->m21 = m21; this->m31 = m31; this->m41 = m41;
    this->m12 = m12; this->m22 = m22; this->m32 = m32; this->m42 = m42;
    this->m13 = m13; this->m23 = m23; this->m33 = m33; this->m43 = m43;
    this->m14 = m14; this->m24 = m24; this->m34 = m34; this->m44 = m44;
}

CMatrix4::CMatrix4(const std::array<float, num_elements>& m)
{
    *this = m;
}

CMatrix4& CMatrix4::operator = (const CMatrix4& m)
{
    std::memcpy(this->m, m.m, sizeof(float) * num_elements);
    return *this;
}

CMatrix4& CMatrix4::operator = (const std::array<float, num_elements>& m)
{
    std::memcpy(this->m, m.data(), sizeof(float) * num_elements);
    return *this;
}

CMatrix4& CMatrix4::operator += (const CMatrix4& m)
{
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] += m.m[i][j];
    return *this;
}

CMatrix4& CMatrix4::operator -= (const CMatrix4& m)
{
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] -= m.m[i][j];
    return *this;
}

CMatrix4& CMatrix4::operator *= (const CMatrix4& m)
{
    *this = *this * m;
    return *this;
}

CMatrix4& CMatrix4::operator *= (float scale)
{
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] *= scale;
    return *this;
}

CMatrix4& CMatrix4::operator /= (float divisor)
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] /= divisor;
    return *this;
}

const CMatrix4 CMatrix4::operator + () const
{
    return *this;
}

const CMatrix4 CMatrix4::operator - () const
{
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = -m[i][j];
    return result;
}

const CMatrix4 CMatrix4::operator + (const CMatrix4& m) const
{
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = this->m[i][j] + m.m[i][j];
    return result;
}

const CMatrix4 CMatrix4::operator - (const CMatrix4& m) const
{
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = this->m[i][j] - m.m[i][j];
    return result;
}

const CMatrix4 CMatrix4::operator * (const CMatrix4& m) const
{
    CMatrix4 result;
    for(auto col = 0; col < order; col++){
        for(auto row = 0; row < order; row++){
            result.m[col][row] = 0.0f;
            for(auto i = 0; i < order; i++){
                result.m[col][row] += (this->m)[i][row] * m.m[col][i];
            }
        }
    }
    return result;
}

const CMatrix4 CMatrix4::operator * (float scale) const
{
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = m[i][j] * scale;
    return result;
}

const CMatrix4 CMatrix4::operator / (float divisor) const
{
    assert(std::fabsf(divisor) > 0.0f);    // division by zero.
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[i][j] = m[i][j] / divisor;
    return result;
}

const CMatrix4 operator * (float scale, const CMatrix4& m)
{
    CMatrix4 result;
    for(auto i = 0; i < CMatrix4::order; i++)
        for(auto j = 0; j < CMatrix4::order; j++)
            result.m[i][j] = scale * m.m[i][j];
    return result;
}

float& CMatrix4::operator () (std::size_t row, std::size_t column)
{
    assert(row < order);
    assert(column < order);
    return m[column][row];
}

const float& CMatrix4::operator () (std::size_t row, std::size_t column) const
{
    assert(row < order);
    assert(column < order);
    return m[column][row];
}

float CMatrix4::Trace() const
{
    return m11 + m22 + m33 + m44;
}

Quaternion CMatrix4::ToRotationQuaternion() const
{
    Quaternion result;
    const float tr = Trace();
    if(tr >= 1.0f){
        // |w| Ç™ç≈ëÂ
        result.w = std::sqrtf(tr) * 0.5f;
        const float rcp_4w = 1.0f / (4.0f * result.w);  // 1/4|w|
        result.x = (m32 - m23) * rcp_4w;    // 4wx / 4|w|
        result.y = (m13 - m31) * rcp_4w;    // 4wy / 4|w|
        result.z = (m21 - m12) * rcp_4w;    // 4wz / 4|w|
    }
    else
    if((m11 > m22) && (m11 > m33)){
        // |x| Ç™ç≈ëÂ
        result.x = std::sqrtf(m11 - m22 - m33 + 1.0f) * 0.5f;
        const float rcp_4x = 1.0f / (4.0f * result.x);  // 1/4|x|
        result.y = (m12 + m21) * rcp_4x;    // 4xy / 4|x|
        result.z = (m13 + m31) * rcp_4x;    // 4xz / 4|x|
        result.w = (m32 - m23) * rcp_4x;    // 4wx / 4|x|
    }
    else
    if((m22 > m33)){
        // |y| Ç™ç≈ëÂ
        result.y = std::sqrtf(m22 - m33 - m11 + 1.0f) * 0.5f;
        const float rcp_4y = 1.0f / (4.0f * result.y);  // 1/4|y|
        result.x = (m12 + m21) * rcp_4y;    // 4xy / 4|y|
        result.z = (m32 + m23) * rcp_4y;    // 4yz / 4|y|
        result.w = (m13 - m31) * rcp_4y;    // 4wy / 4|y|
    }
    else{
        // |z| Ç™ç≈ëÂ
        result.z = std::sqrtf(m33 - m11 - m22 + 1.0f) * 0.5f;
        const float rcp_4z = 1.0f / (4.0f * result.z);  // 1/4|z|
        result.x = (m13 + m31) * rcp_4z;    // 4xz / 4|z|
        result.y = (m32 + m23) * rcp_4z;    // 4yz / 4|z|
        result.w = (m21 - m12) * rcp_4z;    // 4wz / 4|z|
    }
    return result;
}

CMatrix4 CMatrix4::Transpose(const CMatrix4& m)
{
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[j][i] = m.m[i][j];
    return result;
}

CMatrix4 CMatrix4::Inverse(const CMatrix4& m)
{
    assert(0);
    return CMatrix4();
}

CMatrix4 CMatrix4::Translation(float x, float y, float z)
{
    CMatrix4 result = CMatrix4::IDENTITY;
    result.m14 = x;
    result.m24 = y;
    result.m34 = z;
    return result;
}

CMatrix4 CMatrix4::Scaling(float x, float y, float z)
{
    CMatrix4 result = CMatrix4::IDENTITY;
    result.m11 = x;
    result.m22 = y;
    result.m33 = z;
    return result;
}

CMatrix4 CMatrix4::RotationX(float angle)
{
    CMatrix4 result = CMatrix4::IDENTITY;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    result.m22 = c;
    result.m23 =-s;
    result.m32 = s;
    result.m33 = c;
    return result;
}

CMatrix4 CMatrix4::RotationY(float angle)
{
    CMatrix4 result = CMatrix4::IDENTITY;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    result.m11 = c;
    result.m13 = s;
    result.m31 =-s;
    result.m33 = c;
    return result;
}

CMatrix4 CMatrix4::RotationZ(float angle)
{
    CMatrix4 result = CMatrix4::IDENTITY;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    result.m11 = c;
    result.m12 =-s;
    result.m21 = s;
    result.m22 = c;
    return result;
}

CMatrix4 CMatrix4::RotationAxis(Vector3 axis, float angle)
{
    assert(almost_equals(1.0f, axis.Magnitude(), 1));    // axis is not an unit vector.
    const float x = axis.x;
    const float y = axis.y;
    const float z = axis.z;
    const float s = std::sinf(angle);
    const float c = std::cosf(angle);
    const float vers = 1.0f - c;
    return CMatrix4(
        x*x*vers+c,   x*y*vers-z*s, x*z*vers+y*s, 0.0f,
        x*y*vers+z*s, y*y*vers+c,   y*z*vers-x*s, 0.0f,
        x*z*vers-y*s, y*z*vers+x*s, z*z*vers+c,   0.0f,
        0.0f,         0.0f,         0.0f,         1.0f);
}

CMatrix4 CMatrix4::LookAt(Vector3 position, Vector3 target, Vector3 up)
{
    Vector3 zaxis = position - target;
    zaxis.Normalize();
    Vector3 xaxis = Vector3::CrossProduct(up, zaxis);
    xaxis.Normalize();
    Vector3 yaxis = Vector3::CrossProduct(zaxis, xaxis);
    return CMatrix4(
        xaxis.x, xaxis.y, xaxis.z, -Vector3::DotProduct(xaxis, position),
        yaxis.x, yaxis.y, yaxis.z, -Vector3::DotProduct(yaxis, position),
        zaxis.x, zaxis.y, zaxis.z, -Vector3::DotProduct(zaxis, position),
        0.0f,    0.0f,    0.0f,    1.0f);
}

CMatrix4 CMatrix4::Perspective(float fovy, float aspectRatio, float near, float far)
{
    assert(fovy > 0.0f);
    assert(aspectRatio > 0.0f);
    assert(far > near);
    const float cot = 1.0f / std::tanf(fovy * 0.5f);
    const float q = 1.0f / (far - near);
    return CMatrix4(
        cot/aspectRatio, 0.0f, 0.0f,          0.0f,
        0.0f,            cot,  0.0f,          0.0f,
        0.0f,            0.0f, -(far+near)*q, -2.0f*far*near*q,
        0.0f,            0.0f, -1.0f,         0.0f);
}

CMatrix4 CMatrix4::Frustum(float top, float bottom, float left, float right, float near, float far)
{
    assert(top > bottom);
    assert(right > left);
    assert(far > near);
    const float w = 2.0f * near / (right - left);
    const float h = 2.0f * near / (top - bottom);
    const float q = 1.0f / (far - near);
    const float woff = (right + left) / (right - left);
    const float hoff = (top + bottom) / (top - bottom);
    return CMatrix4(
        w,    0.0f, woff,          0.0f,
        0.0f, h,    hoff,          0.0f,
        0.0f, 0.0f, -(far+near)*q, -2.0f*far*near*q,
        0.0f, 0.0f, -1.0f,         0.0f);
}

CMatrix4 CMatrix4::Ortho(float top, float bottom, float left, float right, float near, float far)
{
    assert(top > bottom);
    assert(right > left);
    assert(far > near);
    const float w = 2.0f / (right - left);
    const float h = 2.0f / (top - bottom);
    const float q = 1.0f / (far - near);
    const float woff = -(right + left) / (right - left);
    const float hoff = -(top + bottom) / (top - bottom);
    return CMatrix4(
        w,    0.0f, 0.0f,    woff,
        0.0f, h,    0.0f,    hoff,
        0.0f, 0.0f, -2.0f*q, -(far+near)*q,
        0.0f, 0.0f, 0.0f,    1.0f);
}

std::string CMatrix4::ToString() const
{
    std::string result;
    for(auto i = 0; i < order; i++){
        result.append("CMatrix4[" + std::to_string(i)  + "]{" + std::to_string(m[i][0]));
        for(auto j = 1; j < order; j++){
            result.append(", " + std::to_string(m[i][j]));
        }
        result.append("}\n");
    }
    return result;
}

}}