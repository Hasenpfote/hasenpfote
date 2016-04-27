/*!
* @file cmatrix4.h
* @brief 4x4 column-major matrix class.
* @author Hasenpfote
* @date 2016/04/17
*/
#pragma once
#include <array>

namespace hasenpfote{ namespace math{

class Vector3;
class Quaternion;

class CMatrix4 final
{
public:
    static constexpr auto order = 4;
    static constexpr auto num_elements = order * order;
    union
    {
        struct
        {
            float m11, m21, m31, m41;
            float m12, m22, m32, m42;
            float m13, m23, m33, m43;
            float m14, m24, m34, m44;
        };
        float m[order][order];
    };
public:
    static const CMatrix4 IDENTITY;

public:
/* Constructor */

    CMatrix4() = default;
    CMatrix4(const CMatrix4& m);
    CMatrix4(
        float m11, float m12, float m13, float m14,
        float m21, float m22, float m23, float m24,
        float m31, float m32, float m33, float m34,
        float m41, float m42, float m43, float m44);
    CMatrix4(const std::array<float, num_elements>& m);

/* Destructor */

    ~CMatrix4() = default;

/* Assignment operator */

    CMatrix4& operator = (const CMatrix4& m);
    CMatrix4& operator = (const std::array<float, num_elements>& m);
    CMatrix4& operator += (const CMatrix4& m);
    CMatrix4& operator -= (const CMatrix4& m);
    CMatrix4& operator *= (const CMatrix4& m);
    CMatrix4& operator *= (float scale);
    CMatrix4& operator /= (float divisor);

/* Unary operator */

    const CMatrix4 operator + () const;
    const CMatrix4 operator - () const;

/* Binary operator */

    const CMatrix4 operator + (const CMatrix4& m) const;
    const CMatrix4 operator - (const CMatrix4& m) const;
    const CMatrix4 operator * (const CMatrix4& m) const;
    const CMatrix4 operator * (float scale) const;
    const CMatrix4 operator / (float divisor) const;

    friend const CMatrix4 operator * (float scale, const CMatrix4& m);

/* Subscript operator */
    float& operator () (std::size_t row, std::size_t column);
    const float& operator () (std::size_t row, std::size_t column) const;

/* Operation */

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
     * @return CMatrix4
     */
    static CMatrix4 Transpose(const CMatrix4& m);

    /*!
     * 逆行列を生成.
     * @param[in] m
     * @return CMatrix4
     */
    static CMatrix4 Inverse(const CMatrix4& m);

    /*!
     * 平行移動行列を生成.
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return CMatrix4
     */
    static CMatrix4 Translation(float x, float y, float z);

    /**
     * スケーリング行列を生成.
     * @param[in] x
     * @param[in] y
     * @param[in] z
     * @return CMatrix4
     */
    static CMatrix4 Scaling(float x, float y, float z);

    /*!
     * X 軸周りの回転を表す行列を生成.
     * @param[in] angle an angle in radians.
     * @return CMatrix4
     */
    static CMatrix4 RotationX(float angle);

    /*!
     * Y 軸周りの回転を表す行列を生成.
     * @param[in] angle an angle in radians.
     * @return CMatrix4
     */
    static CMatrix4 RotationY(float angle);

    /*!
     * Z 軸周りの回転を表す行列を生成.
     * @param[in] angle an angle in radians.
     * @return CMatrix4
     */
    static CMatrix4 RotationZ(float angle);

    /*!
     * 任意軸周りの回転を表す行列を生成.
     * @param[in] axis  an unit vector.
     * @param[in] angle an angle in radians.
     * @return CMatrix4
     */
    static CMatrix4 RotationAxis(Vector3 axis, float angle);

    /*!
     * 右手座標系のビュー行列を生成.
     * @param[in] position  位置
     * @param[in] target    注視点
     * @param[in] up        ワールドの上方向
     * @return CMatrix4
     */
    static CMatrix4 LookAt(Vector3 position, Vector3 target, Vector3 up);

    /*!
     * 右手座標系の射影行列を生成(like an OpenGL).
     * @param[in] fovy          total field of view in the YZ plane.(an angle in radians.)
     * @param[in] aspectRatio   aspect ratio of view window.(width:height)
     * @param[in] near          positive distance from camera to near clipping plane.
     * @param[in] far           positive distance from camera to far clipping plane.
     * @return CMatrix4
     */
    static CMatrix4 Perspective(float fovy, float aspectRatio, float near, float far);

    /*!
     * 右手座標系の射影行列を生成(like an OpenGL).
     * @param[in] top       top of view volume at the near clipping plane.
     * @param[in] bottom    bottom of view volume at the near clipping plane.
     * @param[in] left      left of view volume at the near clipping plane.
     * @param[in] right     right of view volume at the near clipping plane.
     * @param[in] near      positive distance from camera to near clipping plane.
     * @param[in] far       positive distance from camera to far clipping plane.
     * @return CMatrix4
     */
    static CMatrix4 Frustum(float top, float bottom, float left, float right, float near, float far);

    /*!
     * 右手座標系の正射影行列を生成(like an OpenGL).
     * @param[in] top       top of parallel view volume.
     * @param[in] bottom    bottom of parallel view volume.
     * @param[in] left      left of parallel view volume.
     * @param[in] right     right of parallel view volume.
     * @param[in] near      positive distance from camera to near clipping plane.
     * @param[in] far       positive distance from camera to far clipping plane.
     * @return CMatrix4
     */
    static CMatrix4 Ortho(float top, float bottom, float left, float right, float near, float far);

/* Debug */
    friend std::ostream& operator<<(std::ostream& os, const CMatrix4& m);
};

/* Inline */

}}