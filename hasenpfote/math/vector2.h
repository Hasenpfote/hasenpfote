/*!
 * @file vector2.h
 * @brief 2d vector class.
 * @author Hasenpfote
 * @date 2016/04/16
 */
#pragma once
#include <array>

namespace hasenpfote{ namespace math{

class Vector2 final
{
public:
    union
    {
        struct
        {
            float x;
            float y;
        };
        std::array<float, 2> v;
    };
public:
    static const Vector2 ZERO;  //!< ゼロベクトル.
    static const Vector2 E_X;   //!< 基底ベクトル.
    static const Vector2 E_Y;   //!< 基底ベクトル.

public:
/* Constructor */

    Vector2() = default;
    Vector2(const Vector2& v);
    Vector2(float x, float y);
    Vector2(const std::array<float, 2>& v);

/* Destructor */

    ~Vector2() = default;

/* Assignment operator */

    Vector2& operator = (const Vector2& v);
    Vector2& operator = (const std::array<float, 2>& v);
    Vector2& operator += (const Vector2& v);
    Vector2& operator -= (const Vector2& v);
    Vector2& operator *= (float scale);
    Vector2& operator /= (float divisor);

/* Unary operator */

    const Vector2 operator + () const;
    const Vector2 operator - () const;

/* Binary operator */

    const Vector2 operator + (const Vector2& v) const;
    const Vector2 operator - (const Vector2& v) const;
    const Vector2 operator * (float scale) const;
    const Vector2 operator / (float divisor) const;

    friend const Vector2 operator * (float scale, const Vector2& v);

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
    Vector2 Normalized() const;

    /*!
     * 符号を反転する.
     */
    void Negate();

/* Static */

    /*!
     * 2 つのベクトルの内積を計算する.
     * @return 内積.
     */
    static float DotProduct(const Vector2& a, const Vector2& b);

    /*!
     * 2 つのベクトルの外積を計算する.
     * @return 外積.
     */
    static float CrossProduct(const Vector2& a, const Vector2& b);

    /*!
     * 2 つのベクトル間の角度を計算する.
     * @return The angle in radians.
     */
    static float Angle(const Vector2& a, const Vector2& b);

    /*!
     * 2 つのベクトル間を線形補間する.
     * @param[in] a 任意のベクトル.
     * @param[in] b 任意のベクトル.
     * @param[in] t [0,1] の補間パラメータ.
     * @return 線形補間されたベクトル.
     */
    static Vector2 Lerp(const Vector2& a, const Vector2& b, float t);

    /*!
     * 2 つのベクトル間を球面線形補間する.
     * @param[in] a 大きさ 1 のベクトル.
     * @param[in] b 大きさ 1 のベクトル.
     * @param[in] t [0,1] の補間パラメータ.
     * @return 球面線形補間されたベクトル.
     */
    static Vector2 Slerp(const Vector2& a, const Vector2& b, float t);

    /*!
     * 重心座標の点を計算する.
     * <p>\f$(1 - f - g)\vec{a} + f\vec{b} + g\vec{c}\f$</p>
     * <p>\f$f\geq 0\f$ かつ \f$g\geq 0\f$ かつ \f$(1 - f - g)\geq 0\f$ ならば三角形内の点.</p>
     * <p>\f$f = 0\f$ かつ \f$g\geq 0\f$ かつ \f$(1 - f - g)\geq 0\f$ ならば線分 \f$\overline{ac}\f$ 上の点.</p>
     * <p>\f$f\geq 0\f$ かつ \f$g = 0\f$ かつ \f$(1 - f - g)\geq 0\f$ ならば線分 \f$\overline{ab}\f$ 上の点.</p>
     * <p>\f$f\geq 0\f$ かつ \f$g\geq 0\f$ かつ \f$(1 - f - g) = 0\f$ ならば線分 \f$\overline{bc}\f$ 上の点.</p>
     * @param[in] a 三角形の頂点 a.
     * @param[in] b 三角形の頂点 b. 
     * @param[in] c 三角形の頂点 c. 
     * @param[in] f 頂点 b に対する重み係数.
     * @param[in] g 頂点 c に対する重み係数.
     * @return 重心座標の点.
     */
    static Vector2 BaryCentric(const Vector2& a, const Vector2& b, const Vector2& c, float f, float g);

    /*!
     * 2 つのベクトルが垂直か.
     * @param[in] a 任意のベクトル.
     * @param[in] b 任意のベクトル.
     * @retval true 垂直である.
     * @retval false 垂直ではない.
     */
    static bool IsPerpendicular(const Vector2& a, const Vector2& b);

    /*!
     * 2 つのベクトルが平行か.
     * @param[in] a 大きさ 1 のベクトル.
     * @param[in] b 大きさ 1 のベクトル.
     * @retval true 平行である.
     * @retval false 平行ではない.
     */
    static bool IsParallel(const Vector2& a, const Vector2& b);

/* Debug */
    std::string ToString() const;
};

/* Inline */

inline Vector2::Vector2(const Vector2& v)
{
    x = v.x;
    y = v.y;
}

inline Vector2::Vector2(float x, float y)
{
    this->x = x;
    this->y = y;
}

inline Vector2::Vector2(const std::array<float, 2>& v)
{
    x = v[0];
    y = v[1];
}

inline Vector2& Vector2::operator = (const Vector2& v)
{
    x = v.x;
    y = v.y;
    return *this;
}

inline Vector2& Vector2::operator = (const std::array<float, 2>& v)
{
    x = v[0];
    y = v[1];
    return *this;
}

inline const Vector2 Vector2::operator + () const
{
    return *this;
}

inline const Vector2 Vector2::operator - () const
{
    return Vector2(-x, -y);
}

}}