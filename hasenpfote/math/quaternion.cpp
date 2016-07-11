#include "../assert.h"
#include <string>
#include "utility.h"
#include "vector3.h"
#include "cmatrix4.h"
#include "rmatrix4.h"
#include "axis_angle.h"
#include "quaternion.h"

namespace hasenpfote{ namespace math{

const Quaternion Quaternion::IDENTITY = Quaternion(1.0f, 0.0f, 0.0f, 0.0f);

Quaternion::Quaternion(const Quaternion& q)
    : w(q.w), x(q.x), y(q.y), z(q.z)
{
}

Quaternion::Quaternion(float w, float x, float y, float z)
    : w(w), x(x), y(y), z(z)
{
}

Quaternion::Quaternion(float s, const Vector3& v)
    : Quaternion(s, v.GetX(), v.GetY(), v.GetZ())
{
}

Quaternion::Quaternion(const Array& q)
    : Quaternion(q[0], q[1], q[2], q[3])
{
}

Quaternion::Quaternion(const AxisAngle& a)
{
    *this = RotationAxis(a.GetAxis(), a.GetAngle());
}

Quaternion& Quaternion::operator = (const Quaternion& q)
{
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;
    return *this;
}

Quaternion& Quaternion::operator = (const Array& q)
{
    w = q[0];
    x = q[1];
    y = q[2];
    z = q[3];
    return *this;
}

Quaternion& Quaternion::operator += (const Quaternion& q)
{
    w += q.w;
    x += q.x;
    y += q.y;
    z += q.z;
    return *this;
}

Quaternion& Quaternion::operator -= (const Quaternion& q)
{
    w -= q.w;
    x -= q.x;
    y -= q.y;
    z -= q.z;
    return *this;
}

Quaternion& Quaternion::operator *= (const Quaternion& q)
{
    const Quaternion temp(*this);
    w = temp.w * q.w - (temp.x * q.x + temp.y * q.y + temp.z * q.z);
    x = temp.w * q.x + q.w * temp.x + (temp.y * q.z - temp.z * q.y);
    y = temp.w * q.y + q.w * temp.y + (temp.z * q.x - temp.x * q.z);
    z = temp.w * q.z + q.w * temp.z + (temp.x * q.y - temp.y * q.x);
    return *this;
}

Quaternion& Quaternion::operator *= (float scale)
{
    w *= scale;
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;
}

Quaternion& Quaternion::operator /= (float divisor)
{
    ASSERT_MSG(std::fabsf(divisor) > 0.0f, "Division by zero.");
    w /= divisor;
    x /= divisor;
    y /= divisor;
    z /= divisor;
    return *this;
}

float Quaternion::NormSquared() const
{
    return w * w + x * x + y * y + z * z;
}

float Quaternion::Norm() const
{
    return std::sqrtf(NormSquared());
}

float Quaternion::NormV() const
{
    return std::sqrtf(x * x + y * y + z * z);
}

void Quaternion::Normalize()
{
    const float n = Norm();
    ASSERT_MSG(n > 0.0f, "Division by zero.");
    w /= n;
    x /= n;
    y /= n;
    z /= n;
}

Quaternion Quaternion::Normalized() const
{
    const float n = Norm();
    ASSERT_MSG(n > 0.0f, "Division by zero.");
    return Quaternion(w / n, x / n, y / n, z / n);
}

CMatrix4 Quaternion::ToRotationCMatrix() const
{
    const float x_sq = x * x;
    const float y_sq = y * y;
    const float z_sq = z * z;
    const float xy = x * y;
    const float yz = y * z;
    const float xz = x * z;
    const float wx = w * x;
    const float wy = w * y;
    const float wz = w * z;
    return CMatrix4(
        1.0f-2.0f*(y_sq+z_sq), 2.0f*(xy-wz),          2.0f*(xz+wy),          0.0f,
        2.0f*(xy+wz),          1.0f-2.0f*(x_sq+z_sq), 2.0f*(yz-wx),          0.0f,
        2.0f*(xz-wy),          2.0f*(yz+wx),          1.0f-2.0f*(x_sq+y_sq), 0.0f,
        0.0f,                  0.0f,                  0.0f,                  1.0f);
}

RMatrix4 Quaternion::ToRotationRMatrix() const
{
    const float x_sq = x * x;
    const float y_sq = y * y;
    const float z_sq = z * z;
    const float xy = x * y;
    const float yz = y * z;
    const float xz = x * z;
    const float wx = w * x;
    const float wy = w * y;
    const float wz = w * z;
    return RMatrix4(
        1.0f-2.0f*(y_sq+z_sq), 2.0f*(xy+wz),          2.0f*(xz-wy),          0.0f,
        2.0f*(xy-wz),          1.0f-2.0f*(x_sq+z_sq), 2.0f*(yz+wx),          0.0f,
        2.0f*(xz+wy),          2.0f*(yz-wx),          1.0f-2.0f*(x_sq+y_sq), 0.0f,
        0.0f,                  0.0f,                  0.0f,                  1.0f);
}

AxisAngle Quaternion::ToAxisAngle() const
{
    ASSERT_MSG(almost_equals(1.0f, Norm(), 1), "Quaternion is not an unit quaternion.");
    AxisAngle result;
    const float i = NormV();
    if(i > 0.0f){    // TODO: 少し余裕を持たせる
        float rcp_i = 1.0f / i;
        result.SetAxis(Vector3(x * rcp_i, y * rcp_i, z * rcp_i));
        result.SetAngle(2.0f * std::atan2f(i, w));
    }
    else{
        result.SetAxis(Vector3::ZERO);
        result.SetAngle(0.0f);
    }
    return result;
}

Vector3 Quaternion::Rotate(const Vector3& v) const
{
    Quaternion p = *this * Quaternion(0.0f, v) * Conjugate(*this);
    return Vector3(p.x, p.y, p.z);
}

Quaternion Quaternion::Inverse(const Quaternion& q)
{
    const float nsq = q.NormSquared();
    ASSERT_MSG(nsq > 0.0f, "Division by zero.");
    return Quaternion(q.w / nsq, -q.x / nsq, -q.y / nsq, -q.z / nsq);
}

Quaternion Quaternion::Conjugate(const Quaternion& q)
{
    return Quaternion(q.w, -q.x, -q.y, -q.z);
}

float Quaternion::DotProduct(const Quaternion& a, const Quaternion& b)
{
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

Quaternion Quaternion::Ln(const Quaternion& q)
{
    const float i = q.NormV();
    const float phi = std::atan2f(i, q.w);
    const float norm = q.Norm();
    float coef;
    if(i > 0.0){    // TODO: 少し余裕を持たせる
        coef = phi / i;
    }
    else{
        coef = rcp_sinc(phi) / norm;
    }
    return Quaternion(std::logf(norm), coef * q.x, coef * q.y, coef * q.z);
}

Quaternion Quaternion::LnU(const Quaternion& q)
{
    const float i = q.NormV();
    const float phi = std::atan2f(i, q.w);
    const float norm = q.Norm();
    ASSERT_MSG(almost_equals(1.0f, norm, 1), "Not an unit quaternion.");

    float _rcp_sinc;
    if(i > 0.0f){    // TODO: 少し余裕を持たせる
        _rcp_sinc = phi / (i / norm);
    }
    else{
        _rcp_sinc = rcp_sinc(phi);
    }
    return Quaternion(0.0f, _rcp_sinc * q.x, _rcp_sinc * q.y, _rcp_sinc * q.z);
}

Quaternion Quaternion::Exp(const Quaternion& q)
{
    const float i = q.NormV();
    float _sinc;
    if(i > 0.0f) {    // TODO: 少し余裕を持たせる
        _sinc = std::sinf(i) / i;
    }
    else{
        _sinc = sinc(i);
    }
    const float exp_w = std::expf(q.w);
    return Quaternion(exp_w * std::cosf(i), exp_w * q.x * _sinc, exp_w * q.y * _sinc, exp_w * q.z * _sinc);
}

Quaternion Quaternion::ExpP(const Quaternion& q)
{
    ASSERT_MSG(almost_equals(0.0f, q.w, 1), "Not a purely imaginary quaternion.");
    const float i = q.NormV();
    float _sinc;
    if(i > 0.0f){    // TODO: 少し余裕を持たせる
        _sinc = std::sinf(i) / i;
    }
    else{
        _sinc = sinc(i);
    }
    return Quaternion(std::cosf(i), q.x * _sinc, q.y * _sinc, q.z * _sinc);
}

Quaternion Quaternion::Pow(const Quaternion& q, float exponent)
{
    return Exp(Ln(q) * exponent);
}

Quaternion Quaternion::PowU(const Quaternion& q, float exponent)
{
    return ExpP(LnU(q) * exponent);
}

Quaternion Quaternion::RotationAxis(const Vector3& axis, float angle)
{
    ASSERT_MSG(almost_equals(1.0f, axis.Magnitude(), 1), "Axis is not an unit quaternion.");
    const float half_angle = angle * 0.5f;
    const float s = std::sinf(half_angle);
    return Quaternion(std::cosf(half_angle), axis * s);
}

Quaternion Quaternion::RotationAxis(const AxisAngle& a)
{
    return RotationAxis(a.GetAxis(), a.GetAngle());
}

Quaternion Quaternion::RotationShortestArc(const Vector3& a, const Vector3& b)
{
    const float d = Vector3::DotProduct(a, b);
    const float s = std::sqrtf((1.0f + d) * 2.0f);
    const Vector3 c = Vector3::CrossProduct(a, b);
    return Quaternion(s * 0.5f, c / s);
}

Quaternion Quaternion::RotationalDifference(const Quaternion& a, const Quaternion& b)
{
    ASSERT_MSG(almost_equals(1.0f, a.Norm(), 1), "Not an unit quaternion.");
    ASSERT_MSG(almost_equals(1.0f, b.Norm(), 1), "Not an unit quaternion.");
    return Quaternion::Conjugate(a) * b;
}

Quaternion Quaternion::Lerp(const Quaternion& a, const Quaternion& b, float t)
{
    ASSERT_MSG(t >= 0.0f && t <= 1.0f, "Not in range.");
    return a + t * (b - a);
}

Quaternion Quaternion::Slerp(const Quaternion& a, const Quaternion& b, float t, bool allowFlip)
{
    ASSERT_MSG(t >= 0.0f && t <= 1.0f, "Not in range.");
    bool flipped = false;           // a または b の反転を表す
    float cos_t = Quaternion::DotProduct(a, b);

    if(allowFlip && (cos_t < 0.0f)){// 最小弧で補間を行う
        flipped = true;
        cos_t = -cos_t;
    }

    float fx, fy;
    if(almost_equals(1.0f, std::fabsf(cos_t), 1)){
        // |cosθ| ≈ 1 → sinθ ≈ 0 の時は線形補間に帰着
        fx = 1.0f - t;
        fy = t;
    }
    else{
        const float theta = std::acosf(cos_t);
        const float cosec = 1.0f / std::sinf(theta);
        fx = std::sinf(theta * (1.0f - t)) * cosec;
        fy = std::sinf(theta * t) * cosec;
    }
    if(flipped){
        fy = -fy;
    }
    return fx * a + fy * b;
}

Quaternion Quaternion::Squad(const Quaternion& p, const Quaternion& q, const Quaternion& a, const Quaternion& b, float t)
{
    ASSERT_MSG(t >= 0.0f && t <= 1.0f, "Not in range.");
    return Slerp(Slerp(p, q, t, false), Slerp(a, b, t, false), 2.0f*t*(1.0f-t), false);
}

Quaternion Quaternion::Spline(const Quaternion& prev, const Quaternion& current, const Quaternion& next)
{
    Quaternion conj_cur = Conjugate(current);
    return current * ExpP(-0.25f * (LnU(conj_cur * prev) + LnU(conj_cur * next)));
}

Quaternion operator + (const Quaternion& q)
{
    return q;
}

Quaternion operator - (const Quaternion& q)
{
    return Quaternion(-q.GetW(), -q.GetX(), -q.GetY(), -q.GetZ());
}

Quaternion operator + (const Quaternion& lhs, const Quaternion& rhs)
{
    return Quaternion(lhs) += rhs;
}

Quaternion operator - (const Quaternion& lhs, const Quaternion& rhs)
{
    return Quaternion(lhs) -= rhs;
}

Quaternion operator * (const Quaternion& lhs, const Quaternion& rhs)
{
    return Quaternion(lhs) *= rhs;
}

Quaternion operator * (const Quaternion& q, float scale)
{
    return Quaternion(q) *= scale;
}

Quaternion operator * (float scale, const Quaternion& q)
{
    return Quaternion(q) *= scale;
}

Quaternion operator / (const Quaternion& q, float divisor)
{
    return Quaternion(q) /= divisor;
}

std::ostream& operator << (std::ostream& os, const Quaternion& q)
{
    const auto flags = os.flags();
    os << "Quaternion{" << q.GetW() << ", " << q.GetX() << ", " << q.GetY() << ", " << q.GetZ() << "}";
    os.flags(flags);
    return os;
}

}}