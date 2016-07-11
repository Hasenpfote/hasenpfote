#include <string>
#include <sstream>
#include "../assert.h"
#include "utility.h"
#include "vector3.h"
#include "vector4.h"
#include "quaternion.h"
#include "cmatrix4.h"

namespace hasenpfote{ namespace math{

const CMatrix4 CMatrix4::ZERO = {
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f
};

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

CMatrix4::CMatrix4(const Vector4& v1, const Vector4& v2, const Vector4& v3, const Vector4& v4)
{
    constexpr std::size_t bytes = sizeof(float) * 4;
    std::memcpy(this->m[0], static_cast<const float*>(v1), bytes);
    std::memcpy(this->m[1], static_cast<const float*>(v2), bytes);
    std::memcpy(this->m[2], static_cast<const float*>(v3), bytes);
    std::memcpy(this->m[3], static_cast<const float*>(v4), bytes);
}

CMatrix4::CMatrix4(const std::array<float, num_elements>& m)
{
    *this = m;
}

void CMatrix4::SetRow(std::int32_t row, const Vector4& v)
{
    ASSERT_MSG((row >= 0) && (row < order), "Row index out of bounds.");
    m[0][row] = v.GetX();
    m[1][row] = v.GetY();
    m[2][row] = v.GetZ();
    m[3][row] = v.GetW();
}

void CMatrix4::SetColumn(std::int32_t column, const Vector4& v)
{
    ASSERT_MSG((column >= 0) && (column < order), "Column index out of bounds.");
    m[column][0] = v.GetX();
    m[column][1] = v.GetY();
    m[column][2] = v.GetZ();
    m[column][3] = v.GetW();
}

Vector4 CMatrix4::GetRow(std::int32_t row) const
{
    ASSERT_MSG((row >= 0) && (row < order), "Row index out of bounds.");
    return Vector4(m[0][row], m[1][row], m[2][row], m[3][row]);
}

Vector4 CMatrix4::GetColumn(std::int32_t column) const
{
    ASSERT_MSG((column >= 0) && (column < order), "Column index out of bounds.");
    return Vector4(m[column][0], m[column][1], m[column][2], m[column][3]);
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
    const CMatrix4 temp(*this);
    for(auto col = 0; col < order; col++){
        for(auto row = 0; row < order; row++){
            this->m[col][row] = 0.0f;
            for(auto i = 0; i < order; i++){
                this->m[col][row] += temp.m[i][row] * m.m[col][i];
            }
        }
    }
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
    ASSERT_MSG(std::fabsf(divisor) > 0.0f, "Division by zero.");
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            this->m[i][j] /= divisor;
    return *this;
}

float& CMatrix4::operator () (std::int32_t row, std::int32_t column)
{
    ASSERT_MSG((row >= 0) && (row < order), "Row index out of bounds.");
    ASSERT_MSG((column >= 0) && (column < order), "Column index out of bounds.");
    return m[column][row];
}

const float& CMatrix4::operator () (std::int32_t row, std::int32_t column) const
{
    ASSERT_MSG((row >= 0) && (row < order), "Row index out of bounds.");
    ASSERT_MSG((column >= 0) && (column < order), "Column index out of bounds.");
    return m[column][row];
}

float CMatrix4::Determinant() const
{
    return m11 * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42))
         + m12 * (m21 * (m34 * m43 - m33 * m44) + m23 * (m31 * m44 - m34 * m41) + m24 * (m33 * m41 - m31 * m43))
         + m13 * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41))
         + m14 * (m21 * (m33 * m42 - m32 * m43) + m22 * (m31 * m43 - m33 * m41) + m23 * (m32 * m41 - m31 * m42));
}

float CMatrix4::Trace() const
{
    return m11 + m22 + m33 + m44;
}

Quaternion CMatrix4::ToRotationQuaternion() const
{
    float w, x, y, z;
    const float tr = Trace();
    if(tr >= 1.0f){
        // |w| が最大
        w = std::sqrtf(tr) * 0.5f;
        const float rcp_4w = 1.0f / (4.0f * w);  // 1/4|w|
        x = (m32 - m23) * rcp_4w;    // 4wx / 4|w|
        y = (m13 - m31) * rcp_4w;    // 4wy / 4|w|
        z = (m21 - m12) * rcp_4w;    // 4wz / 4|w|
    }
    else
    if((m11 > m22) && (m11 > m33)){
        // |x| が最大
        x = std::sqrtf(m11 - m22 - m33 + 1.0f) * 0.5f;
        const float rcp_4x = 1.0f / (4.0f * x);  // 1/4|x|
        y = (m12 + m21) * rcp_4x;    // 4xy / 4|x|
        z = (m13 + m31) * rcp_4x;    // 4xz / 4|x|
        w = (m32 - m23) * rcp_4x;    // 4wx / 4|x|
    }
    else
    if((m22 > m33)){
        // |y| が最大
        y = std::sqrtf(m22 - m33 - m11 + 1.0f) * 0.5f;
        const float rcp_4y = 1.0f / (4.0f * y);  // 1/4|y|
        x = (m12 + m21) * rcp_4y;    // 4xy / 4|y|
        z = (m32 + m23) * rcp_4y;    // 4yz / 4|y|
        w = (m13 - m31) * rcp_4y;    // 4wy / 4|y|
    }
    else{
        // |z| が最大
        z = std::sqrtf(m33 - m11 - m22 + 1.0f) * 0.5f;
        const float rcp_4z = 1.0f / (4.0f * z);  // 1/4|z|
        x = (m13 + m31) * rcp_4z;    // 4xz / 4|z|
        y = (m32 + m23) * rcp_4z;    // 4yz / 4|z|
        w = (m21 - m12) * rcp_4z;    // 4wz / 4|z|
    }
    return Quaternion(w, x, y, z);
}

std::string CMatrix4::ToString() const
{
    std::ostringstream oss;
    for(auto i = 0; i < order; i++){
        oss << "CMatrix4[" << i << "]{" << m[i][0];
        for(auto j = 1; j < order; j++){
            oss << ", " << m[i][j];
        }
        oss << "}" << std::endl;
    }
    return oss.str();
}

CMatrix4 CMatrix4::Transpose(const CMatrix4& m)
{
    CMatrix4 result;
    for(auto i = 0; i < order; i++)
        for(auto j = 0; j < order; j++)
            result.m[j][i] = m.m[i][j];
    return result;
}

CMatrix4 CMatrix4::Inverse(const CMatrix4& m, float* determinant)
{
    const float det = m.Determinant();

    if(determinant)
        *determinant = det;

    if(almost_equals(0.0f, std::fabsf(det), 1))
        return CMatrix4::IDENTITY;

    // 余因子行列を求める.
    CMatrix4 adjugate;
    adjugate.m11 = m.m22 * (m.m33 * m.m44 - m.m34 * m.m43) + m.m23 * (m.m34 * m.m42 - m.m32 * m.m44) + m.m24 * (m.m32 * m.m43 - m.m33 * m.m42);
    adjugate.m21 = m.m21 * (m.m34 * m.m43 - m.m33 * m.m44) + m.m23 * (m.m31 * m.m44 - m.m34 * m.m41) + m.m24 * (m.m33 * m.m41 - m.m31 * m.m43);
    adjugate.m31 = m.m21 * (m.m32 * m.m44 - m.m34 * m.m42) + m.m22 * (m.m34 * m.m41 - m.m31 * m.m44) + m.m24 * (m.m31 * m.m42 - m.m32 * m.m41);
    adjugate.m41 = m.m21 * (m.m33 * m.m42 - m.m32 * m.m43) + m.m22 * (m.m31 * m.m43 - m.m33 * m.m41) + m.m23 * (m.m32 * m.m41 - m.m31 * m.m42);

    adjugate.m12 = m.m12 * (m.m34 * m.m43 - m.m33 * m.m44) + m.m13 * (m.m32 * m.m44 - m.m34 * m.m42) + m.m14 * (m.m33 * m.m42 - m.m32 * m.m43);
    adjugate.m22 = m.m11 * (m.m33 * m.m44 - m.m34 * m.m43) + m.m13 * (m.m34 * m.m41 - m.m31 * m.m44) + m.m14 * (m.m31 * m.m43 - m.m33 * m.m41);
    adjugate.m32 = m.m11 * (m.m34 * m.m42 - m.m32 * m.m44) + m.m12 * (m.m31 * m.m44 - m.m34 * m.m41) + m.m14 * (m.m32 * m.m41 - m.m31 * m.m42);
    adjugate.m42 = m.m11 * (m.m32 * m.m43 - m.m33 * m.m42) + m.m12 * (m.m33 * m.m41 - m.m31 * m.m43) + m.m13 * (m.m31 * m.m42 - m.m32 * m.m41);

    adjugate.m13 = m.m12 * (m.m23 * m.m44 - m.m24 * m.m43) + m.m13 * (m.m24 * m.m42 - m.m22 * m.m44) + m.m14 * (m.m22 * m.m43 - m.m23 * m.m42);
    adjugate.m23 = m.m11 * (m.m24 * m.m43 - m.m23 * m.m44) + m.m13 * (m.m21 * m.m44 - m.m24 * m.m41) + m.m14 * (m.m23 * m.m41 - m.m21 * m.m43);
    adjugate.m33 = m.m11 * (m.m22 * m.m44 - m.m24 * m.m42) + m.m12 * (m.m24 * m.m41 - m.m21 * m.m44) + m.m14 * (m.m21 * m.m42 - m.m22 * m.m41);
    adjugate.m43 = m.m11 * (m.m23 * m.m42 - m.m22 * m.m43) + m.m12 * (m.m21 * m.m43 - m.m23 * m.m41) + m.m13 * (m.m22 * m.m41 - m.m21 * m.m42);

    adjugate.m14 = m.m12 * (m.m24 * m.m33 - m.m23 * m.m34) + m.m13 * (m.m22 * m.m34 - m.m24 * m.m32) + m.m14 * (m.m23 * m.m32 - m.m22 * m.m33);
    adjugate.m24 = m.m11 * (m.m23 * m.m34 - m.m24 * m.m33) + m.m13 * (m.m24 * m.m31 - m.m21 * m.m34) + m.m14 * (m.m21 * m.m33 - m.m23 * m.m31);
    adjugate.m34 = m.m11 * (m.m24 * m.m32 - m.m22 * m.m34) + m.m12 * (m.m21 * m.m34 - m.m24 * m.m31) + m.m14 * (m.m22 * m.m31 - m.m21 * m.m32);
    adjugate.m44 = m.m11 * (m.m22 * m.m33 - m.m23 * m.m32) + m.m12 * (m.m23 * m.m31 - m.m21 * m.m33) + m.m13 * (m.m21 * m.m32 - m.m22 * m.m31);
    // 逆行列を求める.
    return (1.0f / det) * adjugate;
}

CMatrix4 CMatrix4::InverseAffineTransformation(const CMatrix4& m, float* determinant)
{
    const float det = m.m11 * (m.m22 * m.m33 - m.m23 * m.m32)
                    + m.m12 * (m.m23 * m.m31 - m.m21 * m.m33)
                    + m.m13 * (m.m21 * m.m32 - m.m22 * m.m31);

    if(determinant)
        *determinant = det;

    if(almost_equals(0.0f, std::fabsf(det), 1))
        return CMatrix4::IDENTITY;

    CMatrix4 result;

    result.m11 = (m.m22 * m.m33 - m.m23 * m.m32) / det;
    result.m21 = (m.m23 * m.m31 - m.m21 * m.m33) / det;
    result.m31 = (m.m21 * m.m32 - m.m22 * m.m31) / det;
    result.m41 = 0.0f;

    result.m12 = (m.m13 * m.m32 - m.m12 * m.m33) / det;
    result.m22 = (m.m11 * m.m33 - m.m13 * m.m31) / det;
    result.m32 = (m.m12 * m.m31 - m.m11 * m.m32) / det;
    result.m42 = 0.0f;

    result.m13 = (m.m12 * m.m23 - m.m13 * m.m22) / det;
    result.m23 = (m.m13 * m.m21 - m.m11 * m.m23) / det;
    result.m33 = (m.m11 * m.m22 - m.m12 * m.m21) / det;
    result.m43 = 0.0f;

    result.m14 = -(result.m11 * m.m14 + result.m12 * m.m24 + result.m13 * m.m34);
    result.m24 = -(result.m21 * m.m14 + result.m22 * m.m24 + result.m23 * m.m34);
    result.m34 = -(result.m31 * m.m14 + result.m32 * m.m24 + result.m33 * m.m34);
    result.m44 = 1.0f;

    return result;
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
    ASSERT_MSG(almost_equals(1.0f, axis.Magnitude(), 1), "Axis is not an unit vector.");
    const float x = axis.GetX();
    const float y = axis.GetY();
    const float z = axis.GetZ();
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
        xaxis.GetX(), xaxis.GetY(), xaxis.GetZ(), -Vector3::DotProduct(xaxis, position),
        yaxis.GetX(), yaxis.GetY(), yaxis.GetZ(), -Vector3::DotProduct(yaxis, position),
        zaxis.GetX(), zaxis.GetY(), zaxis.GetZ(), -Vector3::DotProduct(zaxis, position),
        0.0f, 0.0f, 0.0f, 1.0f);
}

CMatrix4 CMatrix4::Perspective(float fovy, float aspectRatio, float near, float far)
{
    ASSERT(fovy > 0.0f);
    ASSERT(aspectRatio > 0.0f);
    ASSERT(far > near);
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
    ASSERT(top > bottom);
    ASSERT(right > left);
    ASSERT(far > near);
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
    ASSERT(top > bottom);
    ASSERT(right > left);
    ASSERT(far > near);
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

CMatrix4 operator + (const CMatrix4& m)
{
    return m;
}

CMatrix4 operator - (const CMatrix4& m)
{
    return CMatrix4(CMatrix4::ZERO) -= m;
}

CMatrix4 operator + (const CMatrix4& lhs, const CMatrix4& rhs)
{
    return CMatrix4(lhs) += rhs;
}

CMatrix4 operator - (const CMatrix4& lhs, const CMatrix4& rhs)
{
    return CMatrix4(lhs) -= rhs;
}

CMatrix4 operator * (const CMatrix4& lhs, const CMatrix4& rhs)
{
    return CMatrix4(lhs) *= rhs;
}

CMatrix4 operator * (const CMatrix4& m, float scale)
{
    return CMatrix4(m) *= scale;
}

CMatrix4 operator * (float scale, const CMatrix4& m)
{
    return CMatrix4(m) *= scale;
}

CMatrix4 operator / (const CMatrix4& m, float divisor)
{
    return CMatrix4(m) /= divisor;
}

std::ostream& operator << (std::ostream& os, const CMatrix4& m)
{
    const auto flags = os.flags();
    os << m.ToString() << std::endl;
    os.flags(flags);
    return os;
}

}}