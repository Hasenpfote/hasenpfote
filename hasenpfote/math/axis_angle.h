/*!
* @file axis_angle.h
* @brief Axis-angle class.
* @author Hasenpfote
* @date 2016/04/20
*/
#pragma once
#include <string>

namespace hasenpfote{ namespace math{

class Vector3;

class AxisAngle final
{
public:
    Vector3 axis;   //!< an unit vector.
    float angle;    //!< an angle in radians.

public:
/* Constructor */

    AxisAngle() = default;
    AxisAngle(const AxisAngle& a);
    AxisAngle(const Vector3& axis, float angle);

/* Destructor */

    ~AxisAngle() = default;

/* Assignment operator */

    AxisAngle& operator = (const AxisAngle& a);

/* Debug */
    std::string ToString() const;
};

/* Inline */

inline AxisAngle::AxisAngle(const AxisAngle& a)
{
    axis = a.axis;
    angle = a.angle;
}

inline AxisAngle::AxisAngle(const Vector3& axis, float angle)
{
    this->axis = axis;
    this->angle = angle;
}

inline AxisAngle& AxisAngle::operator = (const AxisAngle& a)
{
    axis = a.axis;
    angle = a.angle;
    return *this;
}

inline std::string AxisAngle::ToString() const
{
    return "AxisAngle{" + axis.ToString() + ", " + std::to_string(angle) + "}";
}

}}