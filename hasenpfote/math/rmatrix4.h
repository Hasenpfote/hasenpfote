/*!
* @file rmatrix4.h
* @brief 4x4 row-major matrix class.
* @author Hasenpfote
* @date 2016/04/19
*/
#pragma once
#include <array>

namespace hasenpfote{ namespace math{

class Vector3;
class Quaternion;

class RMatrix4 final
{
    static constexpr auto order = 4;
    static constexpr auto num_elements = order * order;

public:
/* Constructor */

    RMatrix4() = default;
    RMatrix4(const RMatrix4& m);
    RMatrix4(
        float m11, float m12, float m13, float m14,
        float m21, float m22, float m23, float m24,
        float m31, float m32, float m33, float m34,
        float m41, float m42, float m43, float m44);
    explicit RMatrix4(const std::array<float, num_elements>& m);

/* Destructor */

    ~RMatrix4() = default;

/* Assignment operator */

    RMatrix4& operator = (const RMatrix4& m);
    RMatrix4& operator = (const std::array<float, num_elements>& m);
    RMatrix4& operator += (const RMatrix4& m);
    RMatrix4& operator -= (const RMatrix4& m);
    RMatrix4& operator *= (const RMatrix4& m);
    RMatrix4& operator *= (float scale);
    RMatrix4& operator /= (float divisor);

/* Unary operator */

    const RMatrix4 operator + () const;
    const RMatrix4 operator - () const;

/* Binary operator */

    const RMatrix4 operator + (const RMatrix4& m) const;
    const RMatrix4 operator - (const RMatrix4& m) const;
    const RMatrix4 operator * (const RMatrix4& m) const;
    const RMatrix4 operator * (float scale) const;
    const RMatrix4 operator / (float divisor) const;

    friend const RMatrix4 operator * (float scale, const RMatrix4& m);

/* Subscript operator */
    float& operator () (std::size_t row, std::size_t column);
    const float& operator () (std::size_t row, std::size_t column) const;

/* Operation */

    /*!
     * 行列式.
     * @return
     */
    float Determinant() const;

    /*!
     * トレース.
     * @return
     */
    float Trace() const;

    /*!
     * 回転行列から四元数へ変換する.
     * @return Quaternion
     */
    Quaternion ToRotationQuaternion() const;

/* Static */

    /*!
     * 転置行列を生成.
     * @param[in] m
     * @return RMatrix4
     */
    static RMatrix4 Transpose(const RMatrix4& m);

    /*!
     * 逆行列を生成.
     * @param[in] m
     * @prama[out] determinant
     * @return RMatrix4
     */
    static RMatrix4 Inverse(const RMatrix4& m, float* determinant = nullptr);

    /*!
     * アフィン変換の逆行列を生成.
     * @param[in] m affine transformation matrix
     * @param[out] determinant
     * @return RMatrix4
     */
    static RMatrix4 InverseAffineTransformation(const RMatrix4& m, float* determinant = nullptr);

    /*!
     * 平行移動行列を生成.
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return RMatrix4
     */
    static RMatrix4 Translation(float x, float y, float z);

    /*!
     * スケーリング行列を生成.
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return RMatrix4
     */
    static RMatrix4 Scaling(float x, float y, float z);

    /*!
     * X 軸周りの回転を表す行列を生成.
     * @param[in] angle an angle in radians.
     * @return RMatrix4
     */
    static RMatrix4 RotationX(float angle);

    /*!
     * Y 軸周りの回転を表す行列を生成.
     * @param[in] angle an angle in radians.
     * @return RMatrix4
     */
    static RMatrix4 RotationY(float angle);

    /*!
     * Z 軸周りの回転を表す行列を生成.
     * @param[in] angle an angle in radians.
     * @return RMatrix4
     */
    static RMatrix4 RotationZ(float angle);

    /*!
     * 任意軸周りの回転を表す行列を生成.
     * @param[in] axis  an unit vector.
     * @param[in] angle an angle in radians.
     * @return RMatrix4
     */
    static RMatrix4 RotationAxis(Vector3 axis, float angle);

    /*!
     * 左手座標系のビュー行列を生成.
     * @param[in] position  位置
     * @param[in] target    注視点
     * @param[in] up        ワールドの上方向
     * @return RMatrix4
     */
    static RMatrix4 LookAt(Vector3 position, Vector3 target, Vector3 up);

    /*!
     * 左手座標系の射影行列を生成(like a Direct3D).
     * @param[in] fovy          total field of view in the YZ plane.(an angle in radians.)
     * @param[in] aspectRatio   aspect ratio of view window.(width:height)
     * @param[in] near          positive distance from camera to near clipping plane.
     * @param[in] far           positive distance from camera to far clipping plane.
     * @return RMatrix4
     */
    static RMatrix4 Perspective(float fovy, float aspectRatio, float near, float far);

    /*!
     * 左手座標系の射影行列を生成(like a Direct3D).
     * @param[in] top       top of view volume at the near clipping plane.
     * @param[in] bottom    bottom of view volume at the near clipping plane.
     * @param[in] left      left of view volume at the near clipping plane.
     * @param[in] right     right of view volume at the near clipping plane.
     * @param[in] near      positive distance from camera to near clipping plane.
     * @param[in] far       positive distance from camera to far clipping plane.
     * @return RMatrix4
     */
    static RMatrix4 Frustum(float top, float bottom, float left, float right, float near, float far);

    /*!
     * 左手座標系の正射影行列を生成(like a Direct3D).
     * @param[in] top       top of parallel view volume.
     * @param[in] bottom    bottom of parallel view volume.
     * @param[in] left      left of parallel view volume.
     * @param[in] right     right of parallel view volume.
     * @param[in] near      positive distance from camera to near clipping plane.
     * @param[in] far       positive distance from camera to far clipping plane.
     * @return RMatrix4
     */
    static RMatrix4 Ortho(float top, float bottom, float left, float right, float near, float far);

/* Debug */
    friend std::ostream& operator<<(std::ostream& os, const RMatrix4& m);

public:
    static const RMatrix4 IDENTITY;

public:
    union
    {
        struct
        {
            float m11, m12, m13, m14;
            float m21, m22, m23, m24;
            float m31, m32, m33, m34;
            float m41, m42, m43, m44;
        };
        float m[order][order];
    };
};

/* Inline */

}}