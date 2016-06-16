/*!
* @file vector4.h
* @brief 4d Vector class.
* @author Hasenpfote
* @date 2016/06/16
*/
#pragma once
#include <array>

namespace hasenpfote{ namespace math{

class CMatrix4;
class RMatrix4;

class Vector4 final
{
public:
    union
    {
        struct
        {
            float x;
            float y;
            float z;
            float w;
        };
        std::array<float, 4> v;
    };
public:
    static const Vector4 ZERO;  //!< ゼロベクトル.

public:
/* Constructor */

    Vector4() = default;
    Vector4(const Vector4& v);
    Vector4(float x, float y, float z, float w);
    Vector4(const std::array<float, 4>& v);

/* Destructor */

    ~Vector4() = default;

/* Assignment operator */

    Vector4& operator = (const Vector4& v);
    Vector4& operator = (const std::array<float, 4>& v);
    Vector4& operator += (const Vector4& v);
    Vector4& operator -= (const Vector4& v);
    Vector4& operator *= (float scale);
    Vector4& operator /= (float divisor);

/* Unary operator */

    const Vector4 operator + () const;
    const Vector4 operator - () const;

/* Binary operator */

    const Vector4 operator + (const Vector4& v) const;
    const Vector4 operator - (const Vector4& v) const;
    const Vector4 operator * (float scale) const;
    const Vector4 operator / (float divisor) const;

    friend const Vector4 operator * (float scale, const Vector4& v);
    friend const Vector4 operator * (const CMatrix4& m, const Vector4& v);
    friend const Vector4 operator * (const Vector4& v, const RMatrix4& m);

/* Operation */

    /*!
     * ベクトルの大きさを計算する.
     * @return ベクトルの大きさ.
     */
    float Magnitude() const;

    /*!
     * ベクトルの大きさの二乗を計算する.
     * return ベクトルの大きさの二乗.
     */
    float MagnitudeSquared() const;

    /*!
     * 正規化する.
     */
    void Normalize();

    /*!
     * 正規化したコピーを返す.
     * @return 正規化したベクトル,
     */
    Vector4 Normalized() const;

    /*!
     * 符号を反転する.
     */
    void Negate();

/* Static */

    /*!
     * 2 つのベクトルの内積を計算する.
     * @return 内積.
     */
    static float DotProduct(const Vector4& a, const Vector4& b);

/* Debug */
    friend std::ostream& operator<<(std::ostream& os, const Vector4& v);
};

/* Inline */

}}